/* Description from GLPK:
 *
 * Numbrix is a logic-based number-placement puzzle.[1]
 * The objective is to fill the grid so that each cell contains
 * digits in sequential order taking a horizontal or vertical
 * path; diagonal paths are not allowed. The puzzle setter
 * provides a grid often with the outer most cells completed.
 *
 * Completed Numbrix puzzles are usually a square of numbers
 * in order from 1 to 64 (8x8 grid) or from 1 to 81 (9x9 grid),
 * following a continuous path in sequence.
 *
 * The modern puzzle was invented by Marilyn vos Savant in 2008
 * and published by Parade Magazine under the name "Numbrix",
 * near her weekly Ask Marilyn article.
 *
 *  [1]  http://en.wikipedia.org/wiki/Numbrix
 */

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const int _ = 0;
	vector<vector<int>> given = {{_, _, _, _, _, _, _, _, _},
	                             {_, 11, 12, 15, 18, 21, 62, 61, _},
	                             {_, 6, _, _, _, _, _, 60, _},
	                             {_, 33, _, _, _, _, _, 57, _},
	                             {_, 32, _, _, _, _, _, 56, _},
	                             {_, 37, _, _, _, _, _, 73, _},
	                             {_, 38, _, _, _, _, _, 72, _},
	                             {_, 43, 44, 47, 48, 51, 76, 77, _},
	                             {_, _, _, _, _, _, _, _, _}};

	int n = given.size();

	minimum::linear::IP ip;
	auto x = ip.add_boolean_cube(n, n, n * n);

	// The given numbers are added as hard constraints.
	for (int i : range(n)) {
		for (int j : range(n)) {
			if (given[i][j] > 0) {
				ip.add_constraint(x[i][j][given[i][j] - 1]);
			}
		}
	}

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

	// Every number (except the first) is neighbor to its
	// predecessor.
	for (int k : range(1, n * n)) {
		for (int i : range(n)) {
			for (int j : range(n)) {
				Sum neighborhood = x[i][j][k];

				if (i > 0) {
					neighborhood -= x[i - 1][j][k - 1];
				}
				if (i < n - 1) {
					neighborhood -= x[i + 1][j][k - 1];
				}
				if (j > 0) {
					neighborhood -= x[i][j - 1][k - 1];
				}
				if (j < n - 1) {
					neighborhood -= x[i][j + 1][k - 1];
				}

				ip.add_constraint(neighborhood <= 0);
			}
		}
	}

	IpToSatSolver solver(bind(minisat_solver, false));
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		for (int i : range(n)) {
			for (int j : range(n)) {
				int value = 0;
				for (int k : range(n * n)) {
					value += x[i][j][k].value() * (k + 1);
				}
				cout << setw(3) << right << value;
			}
			cout << "\n";
		}
	}

	return 0;
}

int main() {
	try {
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
