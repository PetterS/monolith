#include <functional>
#include <future>
#include <random>
#include <set>
#include <tuple>
#include <variant>
using namespace std;

#include <catch.hpp>

#include <minimum/core/sqlite.h>
using namespace minimum::core;

TEST_CASE("SqliteStatement_iterate") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (key STRING, value INT);").execute();

	set<pair<string, int>> expected_data{{"a", 1}, {"b", 2}, {"c", 3}};
	auto insert = db.make_statement("INSERT INTO test VALUES(?1, ?2);");
	for (auto& p : expected_data) {
		insert.execute(p.first, p.second);
	}

	auto select1 = db.make_statement<string_view, int32_t>("SELECT * FROM test;");
	set<pair<string, int>> data;
	for (auto& result : select1.execute()) {
		data.emplace(get<0>(result), get<1>(result));
	}
	CAPTURE(to_string(data));
	CHECK(data == expected_data);
}

TEST_CASE("SqliteStatement_raw_iterator") {
	auto db = SqliteDb::sharedInMemory();
	db.make_statement("CREATE TABLE test_raw_iterator (data TEXT);").execute();
	db.make_statement("INSERT INTO test_raw_iterator VALUES(?1);").execute("petter");
	auto select = "SELECT data from test_raw_iterator;";
	CHECK(db.make_statement<string_view>(select).execute().get<0>() == "petter");

	auto stmt = db.make_statement<string_view>(select);
	auto itr = stmt.execute().begin();
	CHECK(get<0>(*itr) == "petter");
	itr++;
	CHECK_THROWS(get<0>(*itr));  // Invalid iterator.

	// The temporary statement is destroyed before begin is
	// called. Since this is a common thing to do, it should
	// work as expected. Everything works because shared_ptr
	// keeps the statement alive.
	auto execution = db.make_statement<string_view>(select).execute();
	auto itr2 = execution.begin();
	CHECK(get<0>(*itr2) == "petter");
	++itr2;
	CHECK_THROWS(get<0>(*itr2));

	// Same here.
	auto itr3 = db.make_statement<string_view>(select).execute().begin();
	CHECK(get<0>(*itr3) == "petter");
	++itr3;
	CHECK_THROWS(get<0>(*itr3));
}

TEST_CASE("SqliteStatement_types") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	db.make_statement("INSERT INTO test VALUES(?1);").execute(1.0);
	auto num = get<0>(*db.make_statement<double>("SELECT * FROM test;").execute().begin());
	CHECK(num == 1.0);
}

TEST_CASE("SqliteStatement_fail_create") {
	auto db = SqliteDb::inMemory();
	CHECK_THROWS_WITH(db.make_statement("FAIL;").execute(), Catch::Contains("syntax error"));
}

TEST_CASE("SqliteStatement_execute_fail") {
	auto db = SqliteDb::inMemory();
	auto stmt = db.make_statement("CREATE TABLE test (data TEXT);");
	stmt.execute();
	CHECK_THROWS_WITH(stmt.execute(), Catch::Contains("table test already exists"));
}

TEST_CASE("SqliteStatement_execute_fail_num_inputs") {
	auto db = SqliteDb::inMemory();
	auto stmt = db.make_statement("CREATE TABLE test (data TEXT);");
	CHECK_THROWS_WITH(stmt.execute(""), Catch::Contains("Wrong number of inputs"));
}

TEST_CASE("SqliteStatement_execute_fail_num_inputs2") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	auto stmt = db.make_statement("INSERT INTO test VALUES(?1);");
	CHECK_THROWS_WITH(stmt.execute(), Catch::Contains("Wrong number of inputs"));
	CHECK_THROWS_WITH(stmt.execute(1, 2), Catch::Contains("Wrong number of inputs"));
}

TEST_CASE("SqliteStatement_execute_fail_num_outputs") {
	auto db = SqliteDb::inMemory();
	CHECK_THROWS_WITH(db.make_statement<int>("CREATE TABLE test (data TEXT);"),
	                  Catch::Contains("Incorrect number of output parameters"));

	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	CHECK_THROWS_WITH(db.make_statement<>("SELECT COUNT(*) FROM test;"),
	                  Catch::Contains("Incorrect number of output parameters"));
}

TEST_CASE("SqliteStatement_execute_fail_not_found") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	auto stmt = db.make_statement<string_view>("SELECT * FROM test WHERE data = ?1;");
	CHECK_NOTHROW(stmt.execute("Petter"));
	// Check that a runtime_error is thrown and not a logic_error.
	CHECK_THROWS_AS(stmt.execute("Petter").get<0>(), std::runtime_error);
}

TEST_CASE("SqliteStatement_get") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	db.make_statement("INSERT INTO test VALUES(NULL);").execute();

	int size = db.make_statement<int>("SELECT COUNT(*) FROM test;").execute().get<0>();
	CHECK(size == 1);
}

TEST_CASE("SqliteStatement_get_all") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (a TEXT, b TEXT);").execute();
	db.make_statement("INSERT INTO test VALUES(?1, ?2);").execute("A", "B");

	auto select = db.make_statement<string_view, string_view>("SELECT * FROM test;");
	auto t = select.execute().get();
	CHECK(get<0>(t) == "A");
	CHECK(get<1>(t) == "B");
}

TEST_CASE("SqliteStatement_possibly_get") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data INTEGER, extra INTEGER);").execute();
	db.make_statement("INSERT INTO test VALUES(?1, ?2);").execute(11, 22);

	auto result = db.make_statement<int>("SELECT data FROM test WHERE data = ?1;")
	                  .execute(1)
	                  .possibly_get<0>();
	CHECK(!result.has_value());
	result = db.make_statement<int>("SELECT data FROM test WHERE data = ?1;")
	             .execute(11)
	             .possibly_get<0>();
	CHECK(*result == 11);

	auto all_results = db.make_statement<int, int>("SELECT data, extra FROM test WHERE data = ?1;")
	                       .execute(1)
	                       .possibly_get();
	CHECK(!all_results.has_value());
	all_results = db.make_statement<int, int>("SELECT data, extra FROM test WHERE data = ?1;")
	                  .execute(11)
	                  .possibly_get();
	CHECK(get<0>(*all_results) == 11);
	CHECK(get<1>(*all_results) == 22);
}

TEST_CASE("SqliteStatement_get_null_string") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	db.make_statement("INSERT INTO test VALUES(NULL);").execute();

	auto sql = "SELECT * FROM test;";
	CHECK(db.make_statement<int>(sql).execute().get<0>() == 0);
	CHECK_THROWS_WITH(db.make_statement<string_view>(sql).execute().get<0>(),
	                  Catch::Contains("string_view can not hold NULL from the database table."));
	CHECK(db.make_statement<optional<string_view>>(sql).execute().get<0>() == nullopt);
	CHECK(db.make_statement<optional<string>>(sql).execute().get<0>() == nullopt);

	db.make_statement("DELETE FROM test;").execute();
	db.make_statement("INSERT INTO test VALUES(?1);").execute("A");
	CHECK(db.make_statement<string_view>(sql).execute().get<0>() == "A");
	CHECK(db.make_statement<optional<string_view>>(sql).execute().get<0>() == "A");
	CHECK(db.make_statement<optional<string>>(sql).execute().get<0>() == "A");
}

TEST_CASE("SqliteStatement_move") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data TEXT);").execute();
	db.make_statement("INSERT INTO test VALUES(?1);").execute("Petter");
	auto stmt = db.make_statement<int>("SELECT COUNT(*) FROM test;");
	CHECK(stmt.execute().get<0>() == 1);
	auto stmt2 = move(stmt);
	CHECK(stmt2.execute().get<0>() == 1);
	stmt = move(stmt2);
	CHECK(stmt.execute().get<0>() == 1);
}

TEST_CASE("SqliteStatement_read_only") {
	auto db = SqliteDb::sharedInMemory();
	const auto cdb = SqliteDb::sharedInMemory();
	db.make_statement("CREATE TABLE read_only (data TEXT);").execute();
	auto insert = "INSERT INTO read_only VALUES(?1);";
	db.make_statement(insert).execute("petter");
	auto select = "SELECT data from read_only;";
	REQUIRE(db.make_statement<string_view>(select).execute().get<0>() == "petter");
	REQUIRE(cdb.make_statement<string_view>(select).execute().get<0>() == "petter");

	// Can not create a statement that writes to the DB.
	CHECK_THROWS(cdb.make_statement(insert));
}

TEST_CASE("SqliteDb_move_assignment") {
	SqliteDb db = SqliteDb::inMemory();
	// Checks for memory leaks in move assignment with sanitizers.
	db = SqliteDb::inMemory();
}

TEST_CASE("SqliteDb_invalid_file") {
	// Checks that no resources are leaked when the SqliteDb
	// throws an exception. (Use e.g. valgrind.)
	CHECK_THROWS(SqliteDb::fromFile("/tmp/"));
}

TEST_CASE("SqliteStatement_string_types") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data STRING);").execute();
	auto insert = db.make_statement("INSERT INTO test VALUES(?1);");
	insert.execute("Test");
	string s = "Petter";
	const string& sref = s;
	insert.execute(sref);
	char* sdata = s.data();
	insert.execute(sdata);
}

TEST_CASE("SqliteDb_custom_function") {
	auto db = SqliteDb::inMemory();
	int offset = 42;
	auto petter = [&offset](int val) {
		if (val == -1) {
			throw std::runtime_error("Invalid Petter");
		}
		return val * val * val + offset;
	};
	db.add_custom_function<int>("petter", petter);

	int val = db.make_statement<int>("SELECT petter(123);").execute().get<0>();
	CHECK(val == 123 * 123 * 123 + offset);
	CHECK_THROWS_WITH(db.make_statement<int>("SELECT petter(-1);").execute(),
	                  Catch::Contains("Invalid Petter"));
}

TEST_CASE("SqliteDb_custom_function_multiple") {
	auto db = SqliteDb::inMemory();
	tuple<int, double, string, string> result;
	auto petter = [&result](int val0, double val1, optional<string_view> val2, string_view val3) {
		get<0>(result) = val0;
		get<1>(result) = val1;
		get<2>(result) = *val2;
		get<3>(result) = val3;
		return 42;
	};
	db.add_custom_function<int, double, optional<string_view>, string_view>("petter", petter);
	int val2 = db.make_statement<int>("SELECT petter(123, 456.789, 'p', 's');").execute().get<0>();
	CHECK(val2 == 42);
	CHECK(get<0>(result) == 123);
	CHECK(get<1>(result) == 456.789);
	CHECK(get<2>(result) == "p");
	CHECK(get<3>(result) == "s");
}

TEST_CASE("SqliteDb_backup") {
	auto db = SqliteDb::inMemory();
	db.make_statement("CREATE TABLE test (data INTEGER);").execute();
	auto insert = db.make_statement("INSERT INTO test VALUES(?1);");
	insert.execute(101);
	insert.execute(102);
	auto db2 = SqliteDb::inMemory();
	auto count_sql = "SELECT count(*) FROM test;";
	REQUIRE_THROWS(db2.make_statement<int>(count_sql));

	db2.backup_from(db);

	CHECK(db2.make_statement<int>(count_sql).execute().get<0>() == 2);
}

TEST_CASE("SqliteDb_backup_fail") {
	auto db = SqliteDb::inMemory();
	CHECK_THROWS(db.backup_from(db));
}

TEST_CASE("SqliteDb_create_table") {
	auto db = SqliteDb::inMemory();
	db.create_table("test", {{"name", "STRING"}, {"age", "INTEGER"}});
	// Recreating with same columns should be a no-op.
	db.create_table("test", {{"name", "STRING"}, {"age", "INTEGER"}});
	db.create_table("test", {{"name", "STRING"}, {"age", "INTEGER"}});

	{
		auto select_age = db.make_statement<int>("SELECT age FROM test WHERE name = ?1;");
		CHECK_NOTHROW(db.make_statement("INSERT INTO test(name, age) VALUES(?1, ?2);")
		                  .execute("Petter", 100));
		CHECK(select_age.execute("Petter").get<0>() == 100);
		// Select statement needs to be destructed before recreating the table.
	}

	db.create_table("test", {{"name", "STRING"}, {"age", "INTEGER"}, {"info", "STRING"}});
	auto select_age = db.make_statement<int>("SELECT age FROM test WHERE name = ?1;");
	auto select_info =
	    db.make_statement<optional<string>>("SELECT info FROM test WHERE name = ?1;");

	CHECK_NOTHROW(db.make_statement("INSERT INTO test(name, age, info) VALUES(?1, ?2, ?3);")
	                  .execute("Gustav", 101, "A"));
	CHECK(select_age.execute("Petter").get<0>() == 100);
	CHECK_FALSE(select_info.execute("Petter").get<0>().has_value());
	CHECK(select_age.execute("Gustav").get<0>() == 101);
	CHECK(select_info.execute("Gustav").get<0>() == "A");

	// Recreating with fewer columns should be OK. Existing data will be kept.
	db.create_table("test", {{"name", "STRING"}, {"age", "INTEGER"}});
	CHECK(select_age.execute("Petter").get<0>() == 100);
	CHECK_FALSE(select_info.execute("Petter").get<0>().has_value());
	CHECK(select_age.execute("Gustav").get<0>() == 101);
	CHECK(select_info.execute("Gustav").get<0>() == "A");

	// Adding and removing columns is currently not supported.
	CHECK_THROWS_WITH(db.create_table("test", {{"name", "STRING"}, {"data", "STRING"}}),
	                  Catch::Contains("no column named age"));
	// Nothing has broken by the failed creation.
	CHECK(select_age.execute("Petter").get<0>() == 100);
	CHECK_FALSE(select_info.execute("Petter").get<0>().has_value());
	CHECK(select_age.execute("Gustav").get<0>() == 101);
	CHECK(select_info.execute("Gustav").get<0>() == "A");
}

TEST_CASE("SqliteDb_move") {
	auto db = SqliteDb::inMemory();
	auto db2 = move(db);
	db = move(db2);
}

TEST_CASE("SqliteDb_many") {
	auto db = SqliteDb::inMemory();
	db.create_table("test", {{"data", "STRING"}});
	mt19937_64 engine;
	uniform_int_distribution<int> dist(1, 100'000);
	auto insert = db.make_statement("INSERT OR IGNORE INTO test VALUES (?1);");

	db.make_statement("BEGIN TRANSACTION;").execute();
	for (int i = 0; i <= 10'000; ++i) {
		insert.execute(to_string(dist(engine)));
	}
	db.make_statement("END TRANSACTION;").execute();

	CHECK(db.make_statement<int>("SELECT COUNT(*) FROM test;").execute().get<0>() > 5000);
}

TEST_CASE("multithreaded_inserts") {
	// Tests that SQLITE_BUSY is handled correctly.
	const int num_iterations = 100;
	auto db = SqliteDb::sharedInMemory();
	db.create_table("threads", {{"id", "INTEGER"}});
	vector<thread> threads;
	for (int i = 0; i < num_iterations; ++i) {
		threads.emplace_back([i]() {
			auto db = SqliteDb::sharedInMemory();
			db.make_statement("INSERT INTO threads(id) VALUES(?);").execute(i);
		});
	}
	for (auto& t : threads) {
		t.join();
	}
	CHECK(db.make_statement<int>("SELECT COUNT(*) FROM threads;").execute().get<0>()
	      == num_iterations);
}
