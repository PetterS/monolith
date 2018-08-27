
#include <cmath>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/testing.h>
using namespace minimum::core;

#define DESCRIPTION(name) if (tests_to_skip.find(name) == tests_to_skip.end())

namespace minimum {
namespace linear {

auto create_soduku_IP(IP& ip, int n) -> decltype(ip.add_boolean_cube(9, 9, 9)) {
	auto P = ip.add_boolean_cube(n * n, n * n, n * n);

	// Exactly one indicator equal to 1.
	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			Sum k_sum;
			for (int k = 0; k < n * n; ++k) {
				k_sum += P[i][j][k];
			}
			ip.add_constraint(k_sum == 1);
		}
	}

	// All rows have every number.
	for (int i = 0; i < n * n; ++i) {
		for (int k = 0; k < n * n; ++k) {
			Sum row_k_sum;
			for (int j = 0; j < n * n; ++j) {
				row_k_sum += P[i][j][k];
			}
			ip.add_constraint(row_k_sum == 1);
		}
	}

	// All columns have every number.
	for (int j = 0; j < n * n; ++j) {
		for (int k = 0; k < n * n; ++k) {
			Sum col_k_sum;
			for (int i = 0; i < n * n; ++i) {
				col_k_sum += P[i][j][k];
			}
			ip.add_constraint(col_k_sum == 1);
		}
	}

	// The n*n subsquares have every number.
	for (int i1 = 0; i1 < n; ++i1) {
		for (int j1 = 0; j1 < n; ++j1) {
			for (int k = 0; k < n * n; ++k) {
				Sum square_k_sum;
				for (int i2 = 0; i2 < n; ++i2) {
					for (int j2 = 0; j2 < n; ++j2) {
						square_k_sum += P[n * i1 + i2][n * j1 + j2][k];
					}
				}
				ip.add_constraint(square_k_sum == 1);
			}
		}
	}
	return P;
}

// Performs simple smoke test for linear programming.
void simple_linear_programming_tests(const Solver& solver, std::set<std::string> tests_to_skip) {
	DESCRIPTION("sudoku") {
		IP ip;
		constexpr int n = 3;
		auto x = create_soduku_IP(ip, n);

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

		for (int i = 0; i < n * n; ++i) {
			for (int j = 0; j < n * n; ++j) {
				if (given[i][j] != '*') {
					int k = given[i][j] - '1';
					ip.add_constraint(x[i][j][k] == 1);
				}
			}
		}

		minimum_core_assert(solver.solutions(&ip)->get());
		minimum_core_assert(ip.is_feasible());
		minimum_core_assert(ip.is_dual_feasible());
	}

	DESCRIPTION("lp") {
		IP ip;
		auto x = ip.add_variable(IP::Real);
		auto y = ip.add_variable(IP::Real);
		ip.add_constraint(x + y >= 3);
		ip.add_constraint(x >= 0);
		ip.add_constraint(y >= 0);
		ip.add_objective(2 * x + 3 * y);

		minimum_core_assert(solver.solutions(&ip)->get());
		minimum_core_assert(std::abs(x.value() - 3) < 1e-6);
		minimum_core_assert(std::abs(y.value()) < 1e-6);
		minimum_core_assert(ip.is_feasible());
		minimum_core_assert(ip.is_dual_feasible());
	}

	DESCRIPTION("lp_unbounded") {
		IP ip;
		auto x = ip.add_variable(IP::Real);
		auto y = ip.add_variable(IP::Real);
		ip.add_constraint(x + y >= 3);
		ip.add_objective(2 * x + 3 * y);

		minimum_core_assert(!solver.solutions(&ip)->get());
	}

	DESCRIPTION("lp_with_duplicates") {
		IP ip;
		auto x = ip.add_variable(IP::Real);
		auto y = ip.add_variable(IP::Real);
		ip.add_constraint(x + y + x + x >= 3);
		ip.add_constraint(x >= 0);
		ip.add_constraint(2 * y == 0);
		ip.add_objective(x);

		minimum_core_assert(solver.solutions(&ip)->get());
		minimum_core_assert(std::abs(x.value() - 1) < 1e-6);
	}
}

// Performs simple smoke tests for integer programming.
void simple_integer_programming_tests(const Solver& solver, std::set<std::string> tests_to_skip) {
	DESCRIPTION("sudoku") {
		IP ip;
		create_soduku_IP(ip);

		minimum_core_assert(solver.solutions(&ip)->get());
		minimum_core_assert(ip.is_feasible_and_integral());
	}

	DESCRIPTION("infeasible") {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y >= 3);
		minimum_core_assert(!solver.solutions(&ip)->get());
	}

	DESCRIPTION("unbounded") {
		IP ip;
		auto x = ip.add_variable(IP::Integer);
		auto y = ip.add_variable(IP::Integer);
		ip.add_constraint(x + y >= 3);
		ip.add_objective(2 * x + 3 * y);

		minimum_core_assert(!solver.solutions(&ip)->get());
	}

	DESCRIPTION("next_solution") {
		IP ip;
		auto x = ip.add_variable(IP::Boolean);
		auto y = ip.add_variable(IP::Boolean);
		auto z = ip.add_variable(IP::Boolean);
		auto w = ip.add_variable(IP::Boolean);
		ip.add_constraint(x + y + z <= 2);
		ip.add_objective(w);

		auto solutions = solver.solutions(&ip);
		int num_solutions = 0;
		while (solutions->get()) {
			num_solutions++;
		}
		minimum_core_assert(num_solutions == 7, "Got ", num_solutions, " solutions.");
	}

	DESCRIPTION("simple_inequality") {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		ip.add_constraint(x + y + z >= 1);

		auto solutions = solver.solutions(&ip);

		int n = 0;
		while (solutions->get()) {
			n++;
		}
		minimum_core_assert(n == 7);
	}

	DESCRIPTION("simple_equality") {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		ip.add_constraint(x + y + z == 1);

		auto solutions = solver.solutions(&ip);

		int n = 0;
		while (solutions->get()) {
			n++;
		}
		minimum_core_assert(n == 3);
	}

	DESCRIPTION("negative_coefficient") {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		ip.add_constraint(x + y + z + w == 1);
		ip.add_constraint(z == w);

		auto solutions = solver.solutions(&ip);

		int n = 0;
		while (solutions->get()) {
			n++;
		}
		minimum_core_assert(n == 2);
	}

	DESCRIPTION("objective") {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		ip.add_constraint(x + y + z + w == 1);
		ip.add_objective(-x - 2 * z + w);

		auto solutions = solver.solutions(&ip);
		minimum_core_assert(solutions->get());
		minimum_core_assert(x.bool_value() == 0);
		minimum_core_assert(y.bool_value() == 0);
		minimum_core_assert(z.bool_value() == 1);
		minimum_core_assert(w.bool_value() == 0);
		minimum_core_assert(!solutions->get());
	}

	DESCRIPTION("square") {
		const int n = 3;

		IP ip;
		auto x = ip.add_boolean_cube(n, n, n * n);
		int s = 15;

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

		int num_solutions = 0;
		auto solutions = solver.solutions(&ip);
		while (solutions->get()) {
			num_solutions++;
		}
		if (n == 3) {
			minimum_core_assert(num_solutions == 8);
		}
	}

	DESCRIPTION("helper_var_1") {
		IP ip;

		ip.add_boolean();
		auto h1 = ip.add_boolean();
		ip.add_boolean();
		auto h2 = ip.add_boolean();
		ip.add_boolean();
		auto h3 = ip.add_boolean();
		auto h4 = ip.add_boolean();

		ip.mark_variable_as_helper(h1);
		ip.mark_variable_as_helper(h2);
		ip.mark_variable_as_helper(h3);
		ip.mark_variable_as_helper(h4);

		// TODO: Glpk crashes without constraints.
		ip.add_constraint(h1 + h2 <= 2);

		int num_solutions = 0;
		auto solutions = solver.solutions(&ip);
		while (solutions->get()) {
			num_solutions++;
		}
		minimum_core_assert(num_solutions == 8);
	}

	DESCRIPTION("helper_var_2") {
		IP ip;
		auto h1 = ip.add_boolean();
		auto h2 = ip.add_boolean();

		ip.mark_variable_as_helper(h1);
		ip.mark_variable_as_helper(h2);

		// TODO: Glpk crashes without constraints.
		ip.add_constraint(h1 + h2 <= 2);

		auto solutions = solver.solutions(&ip);
		minimum_core_assert(solutions->get());
		minimum_core_assert(!solutions->get());
	}

	DESCRIPTION("exists") {
		// http://zverovich.net/2013/02/08/solving-a-puzzle-with-ampl-and-gecode.html
		IP ip;
		auto Digits = {0, 1, 2, 3, 4};
		auto d = ip.add_vector(5, IP::Integer);

		Sum original = 0;
		for (int i : Digits) {
			ip.add_bounds(0, d[i], 9);
			original += d[i] * pow(10, i);
		}

		for (int i : ip.exists(Digits)) {
			Sum sum_with_erased_digit = 0;
			for (int j : Digits) {
				if (j != i) {
					int exp = j < i ? j : j - 1;
					sum_with_erased_digit += d[j] * pow(10, exp);
				}
			}
			ip.add_constraint(original + sum_with_erased_digit == 41751);
		}

		minimum_core_assert(solver.solutions(&ip)->get());
		minimum_core_assert(original.value() == 37956);
	}
}
}  // namespace linear
}  // namespace minimum
