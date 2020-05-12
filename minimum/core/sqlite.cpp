#include <algorithm>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <thread>
#include <utility>

#include <sqlite3.h>

#include <minimum/core/check.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/sqlite.h>
#include <minimum/core/string.h>

using namespace std;

namespace minimum::core {

SqliteDb::SqliteDb(const char* database_file) {
	// This call will only succeed the first time (before sqlite is initialized). That’s OK.
	sqlite3_config(SQLITE_CONFIG_URI, 1);
	int flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;
	auto result = sqlite3_open_v2(database_file, &c_db, flags, nullptr);
	auto guard = make_scope_guard([this]() { sqlite3_close_v2(c_db); });

	if (result != 0) {
		auto msg = to_string(
		    "Could not open database: ", sqlite3_errmsg(c_db), ". File name: ", database_file, ".");
		check(false, msg);
	}

	guard.dismiss();
}

void SqliteDb::create_function(const char* name,
                               int arity,
                               void* data,
                               void (*function)(sqlite3_context*, int n, sqlite3_value**),
                               void (*deleter)(void*)) {
	minimum_core_assert(
	    sqlite3_create_function_v2(c_db, name, arity, 0, data, function, nullptr, nullptr, deleter)
	    == SQLITE_OK);
}

SqliteDb::SqliteDb(SqliteDb&& other) { *this = move(other); }

SqliteDb& SqliteDb::operator=(SqliteDb&& other) noexcept {
	sqlite3_close_v2(c_db);  // nullptr is OK.
	c_db = other.c_db;
	other.c_db = nullptr;
	return *this;
}

SqliteDb::~SqliteDb() noexcept {
	sqlite3_close_v2(c_db);  // nullptr is OK.
}

void SqliteDb::backup_from(const SqliteDb& other) {
	sqlite3_backup* backup = sqlite3_backup_init(db(), "main", other.db(), "main");
	if (backup == nullptr) {
		check(false, "Backup init failed: ", sqlite3_errmsg(db()));
	}
	at_scope_exit(sqlite3_backup_finish(backup));

	int result = -1;
	while ((result = sqlite3_backup_step(backup, 1)) != SQLITE_DONE) {
		// Can print progress here, sleep a while etc....
		if (result == SQLITE_OK) {
			// OK.
		} else if (result == SQLITE_BUSY || result == SQLITE_LOCKED) {
			this_thread::sleep_for(10ms);
		} else {
			check(false, "Backup step error: ", sqlite3_errstr(result));
		}
	}
}

void SqliteDb::create_table(std::string_view name, const std::vector<Column>& columns) {
	set<string> wanted_columns;
	string create_sql = to_string("CREATE TABLE ", name, "(");
	bool first = true;
	for (auto& column : columns) {
		if (!first) {
			create_sql += ", ";
		}
		wanted_columns.emplace(column.name);
		create_sql += column.name;
		create_sql += " ";
		create_sql += column.type_and_constraints;
		first = false;
	}
	create_sql += ");";

	set<string> existing_columns;
	for (auto& column :
	     make_statement<int, string_view>(to_string("PRAGMA table_info(", name, ");")).execute()) {
		existing_columns.emplace(get<1>(column));
	}
	if (existing_columns.empty()) {
		// Table does not exist. Create it.
		make_statement(create_sql).execute();
		return;
	}

	set<string> missing_columns;
	std::set_difference(wanted_columns.begin(),
	                    wanted_columns.end(),
	                    existing_columns.begin(),
	                    existing_columns.end(),
	                    std::inserter(missing_columns, missing_columns.end()));
	if (missing_columns.empty()) {
		// Table exists with all columns. Nothing to do.
		return;
	}

	// Upgrade the table in a transaction.
	make_statement(to_string("BEGIN TRANSACTION;")).execute();
	auto rollback = make_scope_guard([this]() { make_statement("ROLLBACK;").execute(); });

	// No idea what the output column of this statement is supposed to be, but recent version of
	// sqlite started providing it. Since we assert that the expected number of output parameter
	// match, we declare it here.
	make_statement<std::string>(to_string("ALTER TABLE ", name, " RENAME TO create_table_tmp;"))
	    .execute();
	make_statement(create_sql).execute();
	auto existing_columns_str = join(",", existing_columns);
	make_statement(to_string("INSERT INTO ",
	                         name,
	                         "(",
	                         existing_columns_str,
	                         ") SELECT ",
	                         existing_columns_str,
	                         " FROM create_table_tmp;"))
	    .execute();
	make_statement("DROP TABLE create_table_tmp;").execute();

	make_statement("COMMIT;").execute();
	rollback.dismiss();
}

namespace internal {
UntypedSqliteStatement::UntypedSqliteStatement(const SqliteDb& db,
                                               const std::string& sql,
                                               int expected_num_outputs,
                                               ReadState read_state) {
	auto result = sqlite3_prepare_v2(db.db(), sql.c_str(), -1, &stmt, 0);
	auto guard = make_scope_guard([this]() { sqlite3_finalize(stmt); });

	if (result != SQLITE_OK) {
		check(false, "Compiling SQL failed: ", sqlite3_errmsg(db.db()), ".");
	}

	num_input_parameters = sqlite3_bind_parameter_count(stmt);
	num_output_columns = sqlite3_column_count(stmt);

	// Special handling of PRAGMA, which does not seem to work
	// with sqlite3_column_count.
	if (sql.compare(0, 6, "PRAGMA") == 0) {
		// Disable the check.
		num_output_columns = expected_num_outputs;
	}

	check(expected_num_outputs == num_output_columns,
	      "Incorrect number of output parameters when compiling SQL.",
	      " Expected ",
	      expected_num_outputs,
	      " but got ",
	      num_output_columns,
	      ".");
	if (read_state == ReadState::READONLY) {
		check(sqlite3_stmt_readonly(stmt) != 0,
		      "Compiled statement is not read-only for a constant database.");
	}

	guard.dismiss();
}

UntypedSqliteStatement::~UntypedSqliteStatement() {
	sqlite3_finalize(stmt);  // nullptr is OK.
}

std::int32_t UntypedSqliteStatement::get_int32(int pos) { return sqlite3_column_int(stmt, pos); }

std::string_view UntypedSqliteStatement::get_string(int pos) {
	auto ptr = sqlite3_column_text(stmt, pos);
	check(ptr != nullptr, "string_view can not hold NULL from the database table.");
	return reinterpret_cast<const char*>(ptr);
}

std::optional<std::string_view> UntypedSqliteStatement::get_optional_string(int pos) {
	auto ptr = sqlite3_column_text(stmt, pos);
	if (ptr == nullptr) {
		return nullopt;
	}
	return reinterpret_cast<const char*>(ptr);
}
double UntypedSqliteStatement::get_double(int pos) { return sqlite3_column_double(stmt, pos); }

void UntypedSqliteStatement::reset() { sqlite3_reset(stmt); }

void UntypedSqliteStatement::bind(int pos, std::int32_t value) {
	auto result = sqlite3_bind_int(stmt, pos, value);
	minimum_core_assert(result == SQLITE_OK, sqlite3_errstr(result));
}

void UntypedSqliteStatement::bind(int pos, std::string_view value) {
	auto result = sqlite3_bind_text(stmt, pos, value.data(), value.size(), SQLITE_TRANSIENT);
	minimum_core_assert(result == SQLITE_OK, sqlite3_errstr(result));
}

void UntypedSqliteStatement::bind(int pos, double value) {
	auto result = sqlite3_bind_double(stmt, pos, value);
	minimum_core_assert(result == SQLITE_OK, sqlite3_errstr(result));
}

bool UntypedSqliteStatement::execute_statement_once() {
	while (true) {
		auto result = sqlite3_step(stmt);
		if (result == SQLITE_DONE) {
			return false;
		} else if (result == SQLITE_ROW) {
			return true;
		} else if (result == SQLITE_BUSY) {
			this_thread::sleep_for(10ms);
		} else {
			minimum_core_assert(result != SQLITE_MISUSE, sqlite3_errstr(result));
			check(false,
			      "Execution failed: ",
			      sqlite3_errmsg(sqlite3_db_handle(stmt)),
			      " (",
			      sqlite3_errstr(result),
			      ")."
			      "\n");
		}
	}
}

void set_context_result(sqlite3_context* context, int value) { sqlite3_result_int(context, value); }
void set_context_result(sqlite3_context* context, const std::string& value) {
	sqlite3_result_text(context, value.data(), value.size(), SQLITE_TRANSIENT);
}
}  // namespace internal
}  // namespace minimum::core
