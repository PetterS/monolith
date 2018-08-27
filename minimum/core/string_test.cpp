#include <iomanip>
#include <list>
#include <map>

#include <absl/strings/substitute.h>
#include <catch.hpp>

#include <minimum/core/string.h>

using namespace std;
using namespace minimum::core;

TEST_CASE("to_string") { CHECK(to_string("Test", 12, "test") == "Test12test"); }

TEST_CASE("to_string_pair") {
	pair<int, string> p{123, "Test"};
	CHECK(to_string(p) == "(123, Test)");
}

TEST_CASE("to_string_tuple") {
	auto t = make_tuple(123, "Test", 'A');
	CHECK(to_string(t) == "(123, Test, A)");
}

TEST_CASE("to_string_pair_pair") {
	auto p = make_pair(123, "Test");
	auto pp = make_pair(12.1, p);
	CHECK(to_string(pp) == "(12.1, (123, Test))");
}

TEST_CASE("to_string_vector") {
	vector<int> v{1, 2, 3};
	CHECK(to_string(v) == "[1, 2, 3]");
	v.clear();
	CHECK(to_string(v) == "[]");
}

TEST_CASE("to_string_vector_pair") {
	vector<pair<int, string>> v{{1, "P"}, {2, "S"}};
	CHECK(to_string(v) == "[(1, P), (2, S)]");
}

TEST_CASE("to_string_set") {
	set<int> v{1, 2, 3};
	CHECK(to_string(v) == "{1, 2, 3}");
	v.clear();
	CHECK(to_string(v, v, v) == "{}{}{}");
}

TEST_CASE("to_string_map") {
	map<int, string> m;
	m[1] = "P";
	m[0] = "E";
	m[2] = "T";
	CHECK(to_string(m) == "[(0, E), (1, P), (2, T)]");
	CHECK(to_string(m, " ", map<int, int>{}) == "[(0, E), (1, P), (2, T)] []");
	m.clear();
	CHECK(to_string(m) == "[]");
}

TEST_CASE("to_string_array") {
	array<int, 2> arr = {5, 6};
	CHECK(to_string(arr) == "[5, 6]");
}

TEST_CASE("to_string_setprecision") {
	double d = 1.123456;
	CHECK(to_string(setprecision(2), d) == "1.1");
}

TEST_CASE("to_string_setfill_setw") { CHECK(to_string(setfill('P'), setw(3), 1) == "PP1"); }

TEST_CASE("format_string_1") { CHECK(absl::Substitute("Test $0!", 12) == "Test 12!"); }

TEST_CASE("from_string") {
	CHECK(from_string<int>("42") == 42);
	CHECK(from_string("asd", 42) == 42);
	CHECK_THROWS(from_string<int>("abc"));
}

TEST_CASE("join") {
	std::vector<int> v1 = {1, 2, 3};
	CHECK(join('\t', v1) == "1\t2\t3");
	CHECK(join("x", v1) == "1x2x3");

	std::vector<int> v2;
	CHECK(join('\t', v2) == "");
	CHECK(join("\t", v2) == "");

	CHECK(join('\t', {1, 2, 3}) == "1\t2\t3");
	CHECK(join("\t", {1, 2, 3}) == "1\t2\t3");
}

TEST_CASE("to_string_with_separator") {
	CHECK(to_string_with_separator(100) == "100");
	CHECK(to_string_with_separator(1234) == "1,234");
	CHECK(to_string_with_separator(12356) == "12,356");
	CHECK(to_string_with_separator(12345678901LL) == "12,345,678,901");
	CHECK(to_string_with_separator(-12356) == "-12,356");
}

TEST_CASE("split") { CHECK(split("1,2,4", ',') == vector<string>{"1", "2", "4"}); }

TEST_CASE("split some extra") { CHECK(split("1,2,,4,", ',') == vector<string>{"1", "2", "", "4"}); }

TEST_CASE("split empty") {
	CHECK(split("", ',') == vector<string>{});
	// The two below do not really make sense, but were added after the implementation.
	CHECK(split(",", ',') == vector<string>{""});
	CHECK(split(",,,", ',') == vector<string>{"", "", ""});
}

TEST_CASE("strip") {
	auto strip_test = [&](string s) {
		CHECK(strip(s) == s);
		CHECK(strip(" " + s) == s);
		CHECK(strip(s + " ") == s);
		CHECK(strip("    " + s) == s);
		CHECK(strip(s + "  ") == s);
		CHECK(strip(" " + s + " ") == s);
	};
	strip_test("");
	strip_test("s");
	strip_test("123 567");
}

TEST_CASE("remove_spaces") {
	auto remove_spaces_copy = [](string s) {
		remove_spaces(&s);
		return s;
	};
	CHECK(remove_spaces_copy("") == "");
	CHECK(remove_spaces_copy("a") == "a");
	CHECK(remove_spaces_copy("      a  ") == "a");
	CHECK(remove_spaces_copy("123") == "123");
	CHECK(remove_spaces_copy("  123") == "123");
	CHECK(remove_spaces_copy(" 1 2  3   ") == "123");
}

TEST_CASE("ends_with") {
	CHECK(ends_with("Petter", "tter"));
	CHECK(ends_with("Petter", "r"));
	CHECK(ends_with("Petter", ""));
	CHECK_FALSE(ends_with("Petter", "APetter"));
	CHECK_FALSE(ends_with("Petter", "asdf"));
	CHECK_FALSE(ends_with("Petter", "f"));
}
