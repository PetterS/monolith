#include <fstream>

#include <catch.hpp>

#include <minimum/linear/data/util.h>
#include <minimum/linear/projection_solver.h>
#include <minimum/linear/testing.h>
using namespace minimum::linear;

TEST_CASE("sudoku") {
	IP ip;
	int n = 3;
	auto x = create_soduku_IP(ip, n);

	// http://www.mirror.co.uk/news/weird-news/worlds-hardest-sudoku-puzzle-ever-942299
	// clang-format off
	const char* given[] = {"8**" "***" "***",
	                       "**3" "6**" "***",
	                       "*7*" "*9*" "2**",

	                       "*5*" "**7" "***",
	                       "***" "*45" "7**",
	                       "***" "1**" "*3*",

	                       "**1" "***" "*68",
	                       "**8" "5**" "*1*",
	                       "*9*" "***" "4**"};
	// clang-format on
	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			if (given[i][j] != '*') {
				int k = given[i][j] - '1';
				minimum_core_assert(0 <= k && k < n * n);
				ip.add_constraint(x[i][j][k] == 1);
			}
		}
	}

	ProjectionSolver solver;
	solver.tolerance = 1e-3;
	solver.test_interval = 10'000;
	auto solutions = solver.solutions(&ip);
	REQUIRE(solutions->get());
	CHECK(ip.is_feasible(solver.tolerance));
}

void netlib_test(const std::string& name) {
	std::ifstream fin(data::get_directory() + "/netlib/" + name + ".SIF");
	auto ip = read_MPS(fin);

	ProjectionSolver solver;
	solver.test_interval = 100;
	auto solution = solver.solutions(ip.get());
	REQUIRE(solution->get());
	CHECK(ip->is_feasible(solver.tolerance));
}

TEST_CASE("netlib-AFIRO") { netlib_test("AFIRO"); }

TEST_CASE("netlib-SC50A") { netlib_test("SC50A"); }

TEST_CASE("netlib-SC105") { netlib_test("SC105"); }

TEST_CASE("netlib-ADLITTLE") { netlib_test("ADLITTLE"); }

TEST_CASE("netlib-BLEND") { netlib_test("BLEND"); }
