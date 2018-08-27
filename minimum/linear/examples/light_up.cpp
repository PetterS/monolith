//
// A solver for Light Up puzzles.
// http://en.wikipedia.org/wiki/Light_Up_(puzzle)
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const int _ = -2;
	const int X = -1;
	const vector<vector<int>> given = {{X, _, _, X, _, _, _, _, _, X},
	                                   {_, _, _, _, _, _, _, X, _, _},
	                                   {_, 3, _, _, _, _, 0, _, _, _},
	                                   {_, _, 2, _, _, X, _, _, _, 1},
	                                   {_, _, _, 1, 0, X, _, _, _, _},
	                                   {_, _, _, _, 1, X, X, _, _, _},
	                                   {X, _, _, _, 2, _, _, 2, _, _},
	                                   {_, _, _, X, _, _, _, _, X, _},
	                                   {_, _, 1, _, _, _, _, _, _, _},
	                                   {0, _, _, _, _, _, 1, _, _, 0}};

	const int m = given.size();
	const int n = given.back().size();

	IP ip;
	auto x = ip.add_boolean_grid(m, n);

	// Solution must match the given numbers.
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] == X) {
				ip.add_constraint(!x[i][j]);
			} else if (given[i][j] >= 0) {
				Sum neighborhood = 0;
				if (i > 0) {
					neighborhood += x[i - 1][j];
				}
				if (i < n - 1) {
					neighborhood += x[i + 1][j];
				}
				if (j > 0) {
					neighborhood += x[i][j - 1];
				}
				if (j < n - 1) {
					neighborhood += x[i][j + 1];
				}

				ip.add_constraint(!x[i][j]);
				ip.add_constraint(neighborhood == given[i][j]);
			}
		}
	}

	//
	// Every square is illuminated.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] != _) {
				continue;
			}

			Sum vertical_sum = 0;
			for (int i2 = i - 1; i2 >= 0; i2--) {
				if (given[i2][j] != _) {
					break;
				}
				vertical_sum += x[i2][j];
			}
			for (int i2 = i + 1; i2 < m; i2++) {
				if (given[i2][j] != _) {
					break;
				}
				vertical_sum += x[i2][j];
			}

			Sum horizontal_sum = 0;
			for (int j2 = j - 1; j2 >= 0; j2--) {
				if (given[i][j2] != _) {
					break;
				}
				horizontal_sum += x[i][j2];
			}
			for (int j2 = j + 1; j2 < n; j2++) {
				if (given[i][j2] != _) {
					break;
				}
				horizontal_sum += x[i][j2];
			}

			// Lights do not shine on each other. (This will add the same
			// constraint multiple times, but that does not matter.)
			ip.add_constraint(horizontal_sum + x[i][j] <= 1);
			ip.add_constraint(vertical_sum + x[i][j] <= 1);

			// Every square is illuminated.
			ip.add_constraint(horizontal_sum + vertical_sum + x[i][j] >= 1);
		}
	}

	GlpkSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		cout << "\nSolution:\n";
		for (int i : range(n)) {
			for (int j : range(m)) {
				string value = ".";
				if (x[i][j].bool_value()) {
					value = "X";
				} else if (given[i][j] == X) {
					value = "?";
				} else if (given[i][j] >= 0) {
					value = to_string(given[i][j]);
				}
				cout << setw(2) << right << value;
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
