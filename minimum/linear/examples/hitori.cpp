//
// A solver for Hitori puzzles.
// http://en.wikipedia.org/wiki/Hitori
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const vector<vector<int>> given = {{4, 8, 1, 6, 3, 2, 5, 7},
	                                   {3, 6, 7, 2, 1, 6, 5, 4},
	                                   {2, 3, 4, 8, 2, 8, 6, 1},
	                                   {4, 1, 6, 5, 7, 7, 3, 5},
	                                   {7, 2, 3, 1, 8, 5, 1, 2},
	                                   {3, 5, 6, 7, 3, 1, 8, 4},
	                                   {6, 4, 2, 3, 5, 4, 7, 8},
	                                   {8, 7, 1, 4, 2, 3, 5, 6}};

	const int m = given.size();
	const int n = given.back().size();
	IP ip;

	// The color of each square.
	auto x = ip.add_boolean_grid(m, n);

	//
	// No adjacent black squares.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (i > 0) {
				ip.add_constraint(x[i - 1][j] + x[i][j] >= 1);
			}
			if (j > 0) {
				ip.add_constraint(x[i][j - 1] + x[i][j] >= 1);
			}
		}
	}

	//
	// No duplicate numbers in rows and columns.
	//
	for (int i : range(m)) {
		for (int number : range(1, 10)) {
			Sum same = 0;
			int found = 0;
			for (int j : range(n)) {
				if (given[i][j] == number) {
					found++;
					same += x[i][j];
				}
			}
			if (found > 0) {
				ip.add_constraint(same <= 1);
			}
		}
	}
	for (int j : range(n)) {
		for (int number : range(1, 10)) {
			Sum same = 0;
			int found = 0;
			for (int i : range(m)) {
				if (given[i][j] == number) {
					found++;
					same += x[i][j];
				}
			}
			if (found > 0) {
				ip.add_constraint(same <= 1);
			}
		}
	}

	//
	// White region should be connected.
	//
	// ConnectedComponentFlow<> connected(&ip);
	ConnectedComponentNoFlow<> connected(&ip, m * n / 2);
	for (int i : range(m)) {
		for (int j : range(n)) {
			connected.add_node(x[i][j]);
		}
	}
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (i > 0) {
				connected.add_edge(x[i - 1][j], x[i][j]);
			}
			if (j > 0) {
				connected.add_edge(x[i][j - 1], x[i][j]);
			}
		}
	}

	vector<BooleanVariable> all_nodes;
	for (int i : range(m)) {
		for (int j : range(n)) {
			all_nodes.push_back(x[i][j]);
		}
	}

	connected.set_first_node(all_nodes);
	connected.finish();

	IpToSatSolver solver(bind(minisat_solver, false));
	solver.set_silent(true);
	check(solver.solutions(&ip)->get(), "Solving failed.");
	do {
		cout << "\nSolution:\n";
		for (int i : range(m)) {
			for (int j : range(n)) {
				if (x[i][j].bool_value()) {
					cout << given[i][j] << " ";
				} else {
					cout << "  ";
				}
			}
			cout << "\n";
		}
		cout << "\n";
	} while (false && solver.solutions(&ip)->get());

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
