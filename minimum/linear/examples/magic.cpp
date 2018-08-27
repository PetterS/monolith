//
// A solver for Magic squares.
// http://en.wikipedia.org/wiki/Magic_square
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

	const int n = 3;

	IP ip;
	auto x = ip.add_boolean_cube(n, n, n * n);

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

	auto s = ip.add_variable(IP::Real);

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

	GlpkSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	int num_solutions = 0;
	while (solutions->get()) {
		num_solutions++;
		cout << "\nSolution:\n";
		for (int i : range(n)) {
			for (int j : range(n)) {
				cout << setw(3) << right << value[i][j];
			}
			cout << "\n";
		}
	}

	if (n == 3) {
		check(num_solutions == 8, "Wrong number of solutions.");
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
