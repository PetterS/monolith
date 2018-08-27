//
// A solver for Nurikabe puzzles.
// http://en.wikipedia.org/wiki/Nurikabe_(puzzle)
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const int _ = -1;
	/*
	const vector<vector<int>> given = {{2, _, _, _, _, _, _, _, _, 2},
	                                   {_, _, _, _, _, _, 2, _, _, _},
	                                   {_, 2, _, _, 7, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, 3, _, 3, _},
	                                   {_, _, 2, _, _, _, _, 3, _, _},
	                                   {2, _, _, 4, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _},
	                                   {_, 1, _, _, _, _, 2, _, 4, _}};
	//*/

	// #2 from http://fabpedigree.com/nurikabe/ .
	const vector<vector<int>> given = {{5, _, _, _, _, _, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, 5, _, _, _, _, 4},
	                                   {_, _, 3, _, _, _, _, _, _, _, _, _, _, _, _},
	                                   {_, 3, _, _, _, _, _, 4, _, _, _, _, 5, _, _},
	                                   {_, _, _, _, 3, _, _, _, _, 5, _, _, _, 4, _},
	                                   {_, _, 1, _, _, _, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _, _, _, _, _, _},
	                                   {_, _, _, _, 4, _, 4, _, _, _, _, _, _, _, _},
	                                   {_, _, 4, _, _, _, _, _, 1, _, _, _, 1, _, _},
	                                   {_, _, _, _, _, _, _, _, _, _, 3, _, _, _, _},
	                                   {_, _, _, 5, _, _, _, _, _, 3, _, _, _, 4, _},
	                                   {_, 3, _, _, _, _, _, _, _, _, 1, _, _, _, _},
	                                   {_, _, 7, _, _, _, _, 7, _, _, _, _, _, 1, _},
	                                   {_, _, _, _, _, _, _, _, _, _, _, _, _, _, _},
	                                   {4, _, _, _, _, _, _, _, _, _, _, _, _, _, 1}};
	//*/

	const int m = given.size();
	const int n = given.back().size();
	IP ip;
	double start_time = wall_time();

	struct Given {
		int i;
		int j;
		int value;
	};
	vector<Given> givens;

	int num_givens = 0;
	int givens_sum = 0;
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] != _) {
				num_givens++;
				givens.push_back({i, j, given[i][j]});
				givens_sum += given[i][j];
			}
		}
	}

	auto y = ip.add_boolean_cube(givens.size() + 1, m, n);

	//
	// Each squre is a specific numbered component or stream.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			Sum ij_sum = 0;
			for (int g = 0; g < givens.size() + 1; ++g) {
				ij_sum += y[g][i][j];
			}
			ip.add_constraint(ij_sum == 1);
		}
	}

	//
	// Numbered components need to contain their number and
	// be of the right size.
	//
	for (int g = 0; g < givens.size(); ++g) {
		ip.add_bounds(1, y[g][givens[g].i][givens[g].j], 1);

		Sum g_sum = 0;
		for (int i : range(m)) {
			for (int j : range(n)) {
				g_sum += y[g][i][j];
			}
		}
		ip.add_constraint(g_sum == givens[g].value);
	}

	//
	// Numbered components can not be adjacent.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			for (int g1 = 0; g1 < givens.size(); ++g1) {
				for (int g2 = g1 + 1; g2 < givens.size(); ++g2) {
					if (i > 0) {
						ip.add_constraint(y[g1][i - 1][j] + y[g2][i][j] <= 1);
						ip.add_constraint(y[g2][i - 1][j] + y[g1][i][j] <= 1);
					}
					if (j > 0) {
						ip.add_constraint(y[g1][i][j - 1] + y[g2][i][j] <= 1);
						ip.add_constraint(y[g2][i][j - 1] + y[g1][i][j] <= 1);
					}
				}
			}
		}
	}

	//
	// All numbered components are connected.
	//
	for (int g = 0; g < givens.size(); ++g) {
		if (givens[g].value <= 1) {
			continue;
		}
		ConnectedComponentNoFlow<> connected(&ip, givens[g].value);
		for (int i : range(m)) {
			for (int j : range(n)) {
				connected.add_node(y[g][i][j]);
			}
		}

		for (int i : range(m)) {
			for (int j : range(n)) {
				if (i > 0) {
					connected.add_edge(y[g][i - 1][j], y[g][i][j]);
				}
				if (j > 0) {
					connected.add_edge(y[g][i][j - 1], y[g][i][j]);
				}
			}
		}

		connected.finish();
	}

	//
	// The stream is connected.
	//
	int max_stream_length = m * n - givens_sum;
	ConnectedComponentNoFlow<> connected(&ip, max_stream_length);
	for (int i : range(m)) {
		for (int j : range(n)) {
			connected.add_node(y[givens.size()][i][j]);
		}
	}
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (i > 0) {
				connected.add_edge(y[givens.size()][i - 1][j], y[givens.size()][i][j]);
			}
			if (j > 0) {
				connected.add_edge(y[givens.size()][i][j - 1], y[givens.size()][i][j]);
			}
		}
	}
	connected.finish();

	//
	// No 2×2 streams.
	//
	for (int i : range(m - 1)) {
		for (int j : range(n - 1)) {
			ip.add_constraint(y[givens.size()][i][j] + y[givens.size()][i + 1][j]
			                      + y[givens.size()][i][j + 1] + y[givens.size()][i + 1][j + 1]
			                  <= 3);
		}
	}

	//
	// Redundant constraints that speed things up.
	//
	for (int g = 0; g < givens.size(); ++g) {
		auto i0 = givens[g].i;
		auto j0 = givens[g].j;
		for (int i : range(m)) {
			for (int j : range(n)) {
				auto d = abs(i0 - i) + abs(j0 - j);
				if (d >= givens[g].value) {
					ip.add_bounds(0, y[g][i][j], 0);
				}
			}
		}
	}

	double elapsed_time = wall_time() - start_time;
	cout << "Building model: " << elapsed_time << "s.\n";
	start_time = wall_time();

	// Provide starting guesses by solving the LP relaxation.
	IPSolver ip_solver;
	check(ip_solver.solve_relaxation(&ip), "Could not solve relaxation.");

	IpToSatSolver solver(bind(minisat_solver, false));
	check(solver.solutions(&ip)->get(), "Could not solve problem.");

	elapsed_time = wall_time() - start_time;
	cout << "Solving model: " << elapsed_time << "s.\n";

	cout << "\nSolution:\n";
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (y[givens.size()][i][j].bool_value()) {
				cout << "X ";
			} else if (given[i][j] != _) {
				cout << given[i][j] << " ";
			} else {
				cout << "  ";
			}
		}
		cout << "\n";
	}
	cout << "\n";

	return 0;
}

int main(int argc, char* argv[]) {
	try {
		gflags::ParseCommandLineFlags(&argc, &argv, true);
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
