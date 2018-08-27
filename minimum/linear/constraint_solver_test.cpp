#include <iostream>
#include <numeric>

#include <catch.hpp>

#include <minimum/core/range.h>
#include <minimum/linear/constraint_solver.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/testing.h>
using namespace minimum::core;
using namespace minimum::linear;

TEST_CASE("basic") {
	ConstraintSolver solver;
	CHECK_NOTHROW(simple_integer_programming_tests(
	    solver, {"unbounded", "helper_var_1", "helper_var_2", "exists"}));
}

TEST_CASE("square") {
	const int n = 3;

	IP ip;
	auto x = ip.add_boolean_cube(n, n, n * n);
	auto s = ip.add_variable(IP::Integer);
	ip.add_bounds(1, s, 50);

	// Every square contains exactly one number.
	for (int i : range(n)) {
		for (int j : range(n)) {
			Sum k_sum = 0;
			for (int k : range(n * n)) {
				k_sum += x[i][j][k];
			}
			ip.add_constraint(k_sum == 1);
		}
	}

	// Every number occurs exactly once.
	for (int k : range(n * n)) {
		Sum ij_sum = 0;
		for (int i : range(n)) {
			for (int j : range(n)) {
				ij_sum += x[i][j][k];
			}
		}
		ip.add_constraint(ij_sum == 1);
	}

	// Create a helper variable containing the integer at (i, j).
	vector<vector<Sum>> value(n, vector<Sum>(n, 0));
	for (int i : range(n)) {
		for (int j : range(n)) {
			for (int k : range(n * n)) {
				value[i][j] += x[i][j][k] * (k + 1);
			}
		}
	}

	for (int i : range(n)) {
		Sum i_sum = 0;
		for (int j : range(n)) {
			i_sum += value[i][j];
		}
		ip.add_constraint(i_sum == s);
	}

	for (int j : range(n)) {
		Sum j_sum = 0;
		for (int i : range(n)) {
			j_sum += value[i][j];
		}
		ip.add_constraint(j_sum == s);
	}

	Sum diag1 = 0;
	Sum diag2 = 0;
	for (int i : range(n)) {
		diag1 += value[i][i];
		diag2 += value[n - i - 1][i];
	}
	ip.add_constraint(diag1 == s);
	ip.add_constraint(diag2 == s);

	ConstraintSolver solver;
	auto solutions = solver.solutions(&ip);
	int num_solutions = 0;
	while (solutions->get()) {
		CHECK(ip.is_feasible_and_integral());
		num_solutions++;
	}
	if (n == 3) {
		CHECK(num_solutions == 8);
	}
}

TEST_CASE("hexagon") {
	IP ip;

	vector<Sum> x(19, 0);
	vector<Sum> number(19, 0);
	for (int i : range(19)) {
		auto xs = ip.add_boolean_vector(19);
		for (int j : range(19)) {
			x[i] += (j + 1) * xs[j];
			number[j] += xs[j];
		}
		ip.add_constraint(accumulate(xs.begin(), xs.end(), Sum(0)) == 1);
	}

	for (int i : range(19)) {
		ip.add_constraint(number[i] == 1);
	}

	auto rhs = ip.add_variable(IP::Integer);
	ip.add_bounds(1, rhs, 100);

	ip.add_constraint(x[0] + x[1] + x[2] == rhs);
	ip.add_constraint(x[11] + x[12] + x[13] + x[3] == rhs);
	ip.add_constraint(x[10] + x[17] + x[18] + x[14] + x[4] == rhs);
	ip.add_constraint(x[9] + x[16] + x[15] + x[5] == rhs);
	ip.add_constraint(x[8] + x[7] + x[6] == rhs);

	ip.add_constraint(x[10] + x[9] + x[8] == rhs);
	ip.add_constraint(x[11] + x[17] + x[16] + x[7] == rhs);
	ip.add_constraint(x[0] + x[12] + x[18] + x[15] + x[6] == rhs);
	ip.add_constraint(x[1] + x[13] + x[14] + x[5] == rhs);
	ip.add_constraint(x[2] + x[3] + x[4] == rhs);

	ip.add_constraint(x[0] + x[11] + x[10] == rhs);
	ip.add_constraint(x[1] + x[12] + x[17] + x[9] == rhs);
	ip.add_constraint(x[2] + x[13] + x[18] + x[16] + x[8] == rhs);
	ip.add_constraint(x[3] + x[14] + x[15] + x[7] == rhs);
	ip.add_constraint(x[4] + x[5] + x[6] == rhs);

	// Break symmetry.
	ip.add_constraint(x[2] == 3);
	ip.add_constraint(x[3] == 17);
	ip.add_constraint(x[4] == 18);

	ConstraintSolver solver;
	auto solutions = solver.solutions(&ip);
	CHECK(solutions->get());
	CHECK(ip.is_feasible_and_integral());

	CHECK(!solutions->get());
}

TEST_CASE("hexagon-integers") {
	IP ip;
	auto x = ip.add_vector(19, IP::Integer);
	for (int i : range(19)) {
		ip.add_bounds(1, x[i], 19);
	}
	auto rhs = 38;
	ip.add_constraint(x[0] + x[1] + x[2] == rhs);
	ip.add_constraint(x[11] + x[12] + x[13] + x[3] == rhs);
	ip.add_constraint(x[10] + x[17] + x[18] + x[14] + x[4] == rhs);
	ip.add_constraint(x[9] + x[16] + x[15] + x[5] == rhs);
	ip.add_constraint(x[8] + x[7] + x[6] == rhs);

	ip.add_constraint(x[10] + x[9] + x[8] == rhs);
	ip.add_constraint(x[11] + x[17] + x[16] + x[7] == rhs);
	ip.add_constraint(x[0] + x[12] + x[18] + x[15] + x[6] == rhs);
	ip.add_constraint(x[1] + x[13] + x[14] + x[5] == rhs);
	ip.add_constraint(x[2] + x[3] + x[4] == rhs);

	ip.add_constraint(x[0] + x[11] + x[10] == rhs);
	ip.add_constraint(x[1] + x[12] + x[17] + x[9] == rhs);
	ip.add_constraint(x[2] + x[13] + x[18] + x[16] + x[8] == rhs);
	ip.add_constraint(x[3] + x[14] + x[15] + x[7] == rhs);
	ip.add_constraint(x[4] + x[5] + x[6] == rhs);

	ConstraintSolver solver;
	auto solutions = solver.solutions(&ip);
	CHECK(solutions->get());
	CHECK(solutions->get());
	CHECK(solutions->get());
	CHECK(ip.is_feasible_and_integral());
}

TEST_CASE("golomb") {
	IP ip;
	int max_length = 12;  // Maximum length of the ruler.
	int num_vars = max_length + 1;
	auto mark = ip.add_boolean_vector(num_vars);
	// Require a mark at position 0.
	ip.add_bounds(1, mark[0], 1);

	// Either maximize the order or find a ruler with a specific one.
	ip.add_objective(-sum(mark));

	// For each distance, require that it occurs at most once.
	for (int dist = 1; dist < num_vars; ++dist) {
		PseudoBoolean distance_occurrences;
		for (int i : range(num_vars)) {
			for (int j : range(i + 1, num_vars)) {
				if (j - i == dist) {
					distance_occurrences += mark[i] * mark[j];
				}
			}
		}
		ip.add_pseudoboolean_constraint(distance_occurrences <= 1);
	}

	ConstraintSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	int n = 0;
	while (solutions->get()) {
		n++;
		CHECK(sum(mark).value() == 5);
	}
	CHECK(n == 18);
}

TEST_CASE("sudoku-2") {
	IP ip;
	create_soduku_IP(ip, 2);
	ConstraintSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	REQUIRE(solutions->get());
	CHECK(ip.is_feasible_and_integral());
}

TEST_CASE("sudoku-2-assumptions") {
	IP ip;
	auto x = create_soduku_IP(ip, 2);

	ip.add_constraint(x[0][0][3]);
	ip.add_constraint(x[0][1][2]);
	ip.add_constraint(x[1][0][1]);
	ip.add_constraint(x[1][1][0]);

	ConstraintSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	REQUIRE(solutions->get());
	CHECK(ip.is_feasible_and_integral());
}

TEST_CASE("sudoku-2-infeasible") {
	IP ip;
	auto x = create_soduku_IP(ip, 2);

	ip.add_constraint(x[0][0][3]);
	ip.add_constraint(x[0][1][2]);
	ip.add_constraint(x[1][0][1]);
	ip.add_constraint(x[1][1][1]);

	ConstraintSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	REQUIRE(!solutions->get());
}

TEST_CASE("sudoku-3") {
	IP ip;
	auto x = create_soduku_IP(ip, 3);

	const char* given[] = {
	    "1**"
	    "***"
	    "7**",
	    "**7"
	    "1*9"
	    "***",
	    "68*"
	    "*7*"
	    "***",

	    "**1"
	    "*9*"
	    "6**",
	    "***"
	    "3**"
	    "*2*",
	    "*4*"
	    "***"
	    "**3",

	    "**8"
	    "*6*"
	    "1**",
	    "5**"
	    "***"
	    "*4*",
	    "***"
	    "**2"
	    "**5"};
	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 9; ++j) {
			if (given[i][j] != '*') {
				int k = given[i][j] - '1';
				ip.add_constraint(x[i][j][k]);
			}
		}
	}

	ConstraintSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	// REQUIRE(solutions->get());
}
