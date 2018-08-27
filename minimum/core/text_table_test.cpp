#include <sstream>
#include <string>

#include <catch.hpp>

#include <minimum/core/text_table.h>

using namespace std;

TEST_CASE("table") {
	// clang-format off
	vector<vector<int>> data = {
		{ 1, 1, 1 },
		{ 1, 2, 2 },
		{ 1, 1, 1 },
	};
	string expected =
		"+--------+\n"
		"| 1  1  1|\n"
		"|  +-----+\n"
		"| 1| 2  2|\n"
		"|  +-----+\n"
		"| 1  1  1|\n"
		"+--------+\n";
	// clang-format on

	ostringstream out;
	print_table(out, data);
	CHECK(out.str() == expected);
}

TEST_CASE("non_square_table") {
	// clang-format off
	vector<vector<int>> data = {
		{ 1, 2 },
		{ 2, 1 },
		{ 2, 1 },
		{ 1, 1 },
	};
	string expected =
		"+--+--+\n"
		"| 1| 2|\n"
		"+--+--+\n"
		"| 2| 1|\n"
		"|  |  |\n"
		"| 2| 1|\n"
		"+--+  |\n"
		"| 1  1|\n"
		"+-----+\n";
	// clang-format on

	ostringstream out;
	print_table(out, data);
	CHECK(out.str() == expected);
}
