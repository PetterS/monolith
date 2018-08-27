#pragma once
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>

#include <sqlite3.h>

#include <minimum/core/check.h>
#include <minimum/core/export.h>

namespace minimum::core {

// There are two main classes in the API:
// Holds the database connection.
class SqliteDb;
// A compiled SQL statement. Created by SqliteDb.
template <typename... OutputTypes>
class SqliteStatement;

namespace internal {
class UntypedSqliteStatement;

// Never created (or even stored) explicity. Holds the result of an
// executed statement.
//
// Either call get<i>(), begin(), or just iterate over the results
// in a range-based for loop.
template <typename... RowTypes>
class SqliteResult {
   public:
	friend class ::minimum::core::internal::UntypedSqliteStatement;
	class Iterator {
	   public:
		friend class SqliteResult;

		const std::tuple<RowTypes...>& operator*() const {
			check(stmt != nullptr, "No result is available.");
			return result;
		}
		bool operator==(const Iterator& rhs) const { return stmt == rhs.stmt; }
		bool operator!=(const Iterator& rhs) const { return !(*this == rhs); }
		void operator++();
		void operator++(int) { ++(*this); }
		operator bool() const { return stmt != nullptr; }

	   private:
		Iterator() {}
		Iterator(std::shared_ptr<UntypedSqliteStatement> stmt_) : stmt(std::move(stmt_)) {}

		template <int start>
		void exctract_all();

		std::shared_ptr<UntypedSqliteStatement> stmt = nullptr;
		std::tuple<RowTypes...> result;
	};
	Iterator begin() { return itr; }
	Iterator end() { return {}; }

	template <int pos>
	auto get() {
		return std::get<pos>(*itr);
	}

	auto get() { return *itr; }

	template <int pos>
	std::optional<typename std::tuple_element<pos, std::tuple<RowTypes...>>::type> possibly_get() {
		if (!itr) {
			return std::nullopt;
		}
		return std::make_optional(std::get<pos>(*itr));
	}

	std::optional<std::tuple<RowTypes...>> possibly_get() {
		if (!itr) {
			return std::nullopt;
		}
		return std::make_optional(*itr);
	}

   private:
	SqliteResult(std::shared_ptr<UntypedSqliteStatement> stmt) : itr(std::move(stmt)) { ++itr; }
	Iterator itr;
};

// Used internally to hold a compiled statement.
class MINIMUM_CORE_API UntypedSqliteStatement
    : public std::enable_shared_from_this<UntypedSqliteStatement> {
   public:
	friend class SqliteDb;
	template <typename... RowTypes>
	friend class SqliteResult;

	// Executes the statement for the specified output types and inputs.
	//
	// Output types are specified explicity, like:
	//
	//    statement.execute<int, int>("test");  // Output is (int, int).
	//
	// Throws an exception if there is an execution error, otherwise
	// returns the result in an SqliteResult.
	template <typename... OutputTypes, typename... InputTypes>
	SqliteResult<OutputTypes...> execute(InputTypes... input) {
		minimum_core_assert(sizeof...(OutputTypes) == num_output_columns);
		check(sizeof...(InputTypes) == num_input_parameters,
		      "Wrong number of inputs in execute.",
		      " Expected ",
		      num_input_parameters,
		      " but got ",
		      sizeof...(InputTypes),
		      ".");
		reset();
		bind_all(input...);
		return SqliteResult<OutputTypes...>(shared_from_this());
	}

	// Not called except internally. Clang has problems making this
	// private.
	bool execute_statement_once();

	~UntypedSqliteStatement();

	enum class ReadState { READONLY, READWRITE };
	UntypedSqliteStatement(const SqliteDb& db,
	                       const std::string& sql,
	                       int expected_num_outputs,
	                       ReadState read_state);

   private:
	UntypedSqliteStatement(const UntypedSqliteStatement&) = delete;

	void reset();

	void bind(int pos, std::int32_t value);
	void bind(int pos, std::string_view value);
	void bind(int pos, double value);

	template <typename T, typename... Ts>
	void bind_all(T t, Ts... ts) {
		int pos = num_input_parameters - sizeof...(Ts);
		bind(pos, t);
		bind_all(ts...);
	}
	void bind_all() {}

	std::int32_t get_int32(int pos);
	std::string_view get_string(int pos);
	std::optional<std::string_view> get_optional_string(int pos);
	double get_double(int pos);

	sqlite3_stmt* stmt = nullptr;
	int num_input_parameters = -1;
	int num_output_columns = -1;
};
}  // namespace internal

// A compiled SQL statement.
//
// Created by SqliteDb.
template <typename... OutputTypes>
class SqliteStatement {
   public:
	SqliteStatement(std::shared_ptr<internal::UntypedSqliteStatement> stmt_)
	    : stmt(std::move(stmt_)) {}

	// Executes the statement for the specified output types and inputs.
	//
	//    statement.execute("test");
	//
	// Throws an exception if there is an execution error, otherwise
	// returns the result in an SqliteResult.
	template <typename... InputTypes>
	internal::SqliteResult<OutputTypes...> execute(InputTypes... input) {
		return stmt->execute<OutputTypes...>(input...);
	}

	SqliteStatement(const SqliteStatement&) = delete;
	SqliteStatement& operator=(const SqliteStatement&) = delete;
	SqliteStatement(SqliteStatement&&) = default;
	SqliteStatement& operator=(SqliteStatement&&) = default;

   private:
	std::shared_ptr<internal::UntypedSqliteStatement> stmt;
};

namespace internal {
// More overloads needed.
MINIMUM_CORE_API void set_context_result(struct sqlite3_context* context, int value);
MINIMUM_CORE_API void set_context_result(struct sqlite3_context* context, const std::string& value);

template <int start, typename... ArgTypes>
void get_arguments(std::tuple<ArgTypes...>& args, sqlite3_value** values);
}  // namespace internal

class MINIMUM_CORE_API SqliteDb {
   public:
	// Creates an sqlite database in a file.
	static SqliteDb fromFile(const std::string& file_name) { return SqliteDb(file_name.c_str()); }
	// Creates a database in memory that is not shared between instances.
	static SqliteDb inMemory() { return SqliteDb(":memory:"); }
	// Creates a database in memory that is shared between instances.
	static SqliteDb sharedInMemory() { return SqliteDb("file::memory:?cache=shared"); }

	SqliteDb(SqliteDb&& other);
	SqliteDb& operator=(SqliteDb&& other);
	~SqliteDb();

	// Compiles an SQL statement. The only way to interact with the database.
	template <typename... OutputTypes>
	SqliteStatement<OutputTypes...> make_statement(const std::string& sql) {
		// Can not use make_shared (private constructor).
		auto read_state = internal::UntypedSqliteStatement::ReadState::READWRITE;
		auto untyped =
		    new internal::UntypedSqliteStatement(*this, sql, sizeof...(OutputTypes), read_state);
		return {std::shared_ptr<internal::UntypedSqliteStatement>(untyped)};
	}

	// Creates a statement that can only read the database. Throws an exception
	// if the statement is not read-only after compilation.
	template <typename... OutputTypes>
	SqliteStatement<OutputTypes...> make_statement(const std::string& sql) const {
		// Can not use make_shared (private constructor).
		auto read_state = internal::UntypedSqliteStatement::ReadState::READONLY;
		auto untyped =
		    new internal::UntypedSqliteStatement(*this, sql, sizeof...(OutputTypes), read_state);
		return {std::shared_ptr<internal::UntypedSqliteStatement>(untyped)};
	}

	struct Column {
		std::string_view name;
		std::string_view type_and_constraints;
	};

	// Creates a new table with the provided columns.
	//  • Table does not exist: it is created.
	//  • Table exists and has all column names: no change.
	//  • Table exists but one column name is missing: table is rebuilt
	//    and data in the existing column names is kept.
	// Does not handle foreign keys, indices etc.
	void create_table(std::string_view name, const std::vector<Column>& columns);

	// Assigns the content of other into this. The other DB may be used
	// continuously throughout the backup process.
	void backup_from(const SqliteDb& other);

	// Adds a custom  function that can be used in subsequent statements.
	// The function should return string or int. Exceptions will be caught.
	template <typename... ArgTypes, typename Function>
	void add_custom_function(const char* name, Function&& function_) {
		constexpr int arity = sizeof...(ArgTypes);
		static_assert(arity > 0, "Function arguments must be manually specified.");

		auto function = new auto(std::forward<Function>(function_));
		auto deleter = [](void* ptr) -> void { delete static_cast<decltype(function)>(ptr); };
		auto wrapper = [](sqlite3_context* context, int n, sqlite3_value** values) -> void {
			constexpr int arity = sizeof...(ArgTypes);
			auto& f = *static_cast<decltype(function)>(sqlite3_user_data(context));
			std::tuple<ArgTypes...> args;
			internal::get_arguments<arity - 1>(args, values);
			try {
				auto result = std::apply(f, args);
				internal::set_context_result(context, result);
			} catch (std::exception& err) {
				sqlite3_result_error(context, err.what(), -1);
			} catch (...) {
				sqlite3_result_error(context, "Unknown exception.", -1);
			}
		};
		create_function(name, arity, function, wrapper, deleter);
	}

	// Access to sqlite internals.
	sqlite3* db() const { return c_db; }

   private:
	SqliteDb(const char* database_file);
	SqliteDb(const SqliteDb&) = delete;

	void create_function(const char* name,
	                     int arity,
	                     void* data,
	                     void (*function)(sqlite3_context*, int n, sqlite3_value**),
	                     void (*deleter)(void*));

	sqlite3* c_db;
};

namespace internal {
template <typename... RowTypes>
void SqliteResult<RowTypes...>::Iterator::operator++() {
	if (stmt == nullptr) {
		// This iterator is already the end iterator.
		return;
	}
	if (!stmt->execute_statement_once()) {
		// Execution succeeded, but there were no results.
		// Iterator becomes the end iterator.
		stmt = nullptr;
		return;
	}
	// There were results; extract them to the tuple.
	exctract_all<int(sizeof...(RowTypes)) - 1>();
}

template <typename... RowTypes>
template <int start>
void SqliteResult<RowTypes...>::Iterator::exctract_all() {
	if constexpr (start <= -1) {
		return;
	} else {
		using T = typename std::tuple_element_t<start, std::tuple<RowTypes...>>;

		constexpr bool is_string =
		    (std::is_same_v<T, std::string_view>) || (std::is_same_v<T, std::string>)
		    || (std::is_same_v<T, const char*>);
		constexpr bool is_optional_string =
		    std::is_same_v<
		        T,
		        std::optional<std::string_view>> || std::is_same_v<T, std::optional<std::string>>;
		constexpr bool is_int32 = std::is_same_v<T, std::int32_t>;

		static_assert(is_string || is_optional_string || is_int32 || std::is_floating_point_v<T>);
		if constexpr (is_string) {
			std::get<start>(result) = stmt->get_string(start);
		} else if constexpr (is_optional_string) {
			std::get<start>(result) = stmt->get_optional_string(start);
		} else if constexpr (is_int32) {
			std::get<start>(result) = stmt->get_int32(start);
		} else if (std::is_floating_point_v<T>) {
			std::get<start>(result) = stmt->get_double(start);
		}
		exctract_all<start - 1>();
	}
}

template <int start, typename... ArgTypes>
void get_arguments(std::tuple<ArgTypes...>& args, sqlite3_value** values) {
	if constexpr (start <= -1) {
		return;
	} else {
		using T = typename std::tuple_element_t<start, std::tuple<ArgTypes...>>;

		constexpr bool is_string =
		    (std::is_same_v<T, std::string_view>) || (std::is_same_v<T, std::string>)
		    || (std::is_same_v<T, const char*>);
		constexpr bool is_optional_string =
		    std::is_same_v<
		        T,
		        std::optional<std::string_view>> || std::is_same_v<T, std::optional<std::string>>;
		constexpr bool is_int32 = std::is_same_v<T, std::int32_t>;

		static_assert(is_string || is_optional_string || is_int32 || std::is_floating_point_v<T>);
		if constexpr (is_string) {
			auto ptr = sqlite3_value_text(values[start]);
			check(ptr != nullptr, "string_view can not hold NULL from the database table.");
			std::get<start>(args) = reinterpret_cast<const char*>(ptr);
		} else if constexpr (is_optional_string) {
			auto ptr = sqlite3_value_text(values[start]);
			if (ptr == nullptr) {
				std::get<start>(args) = std::nullopt;
			} else {
				std::get<start>(args) = reinterpret_cast<const char*>(ptr);
			}
		} else if constexpr (is_int32) {
			std::get<start>(args) = sqlite3_value_int(values[start]);
		} else if constexpr (std::is_floating_point_v<T>) {
			std::get<start>(args) = sqlite3_value_double(values[start]);
		}
		get_arguments<start - 1>(args, values);
	}
}
}  // namespace internal
}  // namespace minimum::core
