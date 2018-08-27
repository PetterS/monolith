// Petter Strandmark 2013
// petter.strandmark@gmail.com

#include <iostream>
#include <stdexcept>
#include <vector>
using std::size_t;
using std::vector;

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
using namespace minimum::core;
using namespace minimum::linear;

//
// This class adds connectivity constraints to the IP class.
// It does so in two different ways, but they both result in
// problems that are quite difficult to solve.
//
class MyIP : public IP {
   public:
	// The first way is based on following the perimeter around the
	// area and counting left and right turns.
	void add_connected_constraint_perimeter(const vector<vector<BooleanVariable>>& grid) {
		auto m = grid.size();
		if (m == 0) return;
		auto n = grid[0].size();

		// Pad grid with zeroes to be able to add constraints everywhere.
		vector<vector<BooleanVariable>> padded_grid;
		for (int i = 0; i < m + 2; ++i) {
			padded_grid.push_back(vector<BooleanVariable>());
			for (int j = 0; j < n + 2; ++j) {
				if (i == 0 || i == m + 1 || j == 0 || j == n + 1) {
					auto zero_var = add_boolean();
					add_bounds(0, zero_var, 0);
					padded_grid.back().push_back(zero_var);
				} else {
					padded_grid.back().push_back(grid[i - 1][j - 1]);
				}
			}
		}

		Sum turn_sum;

		auto add_triple_constraints =
		    [&](const BooleanVariable& x1, const BooleanVariable& x2, const BooleanVariable& y) {
			    auto v1 = add_boolean();
			    // v1  <=>  y && !x1 && !x2.
			    // <=
			    add_constraint(y + (1 - x1) + (1 - x2) <= 2 + v1);
			    // =>
			    add_constraint(implication(v1, !x1));
			    add_constraint(implication(v1, !x2));
			    add_constraint(implication(v1, y));
			    turn_sum += v1;

			    auto v2 = add_boolean();
			    // v2  <=>  !y && x1 && x2.
			    // <=
			    add_constraint((1 - y) + x1 + x2 <= 2 + v2);
			    // =>
			    add_constraint(implication(v2, x1));
			    add_constraint(implication(v2, x2));
			    add_constraint(implication(v2, !y));
			    turn_sum -= v2;
		    };

		for (size_t i = 0; i < m + 1; ++i) {
			for (size_t j = 0; j < n + 1; ++j) {
				//  X
				//  YX
				{
					auto& x1 = padded_grid[i][j];
					auto& x2 = padded_grid[i + 1][j + 1];
					auto& y = padded_grid[i][j + 1];
					add_triple_constraints(x1, x2, y);
				}

				// XY
				//  X
				{
					auto& x1 = padded_grid[i][j];
					auto& x2 = padded_grid[i + 1][j + 1];
					auto& y = padded_grid[i + 1][j];
					add_triple_constraints(x1, x2, y);
				}

				// YX
				// X
				{
					auto& x1 = padded_grid[i + 1][j];
					auto& x2 = padded_grid[i][j + 1];
					auto& y = padded_grid[i][j];
					add_triple_constraints(x1, x2, y);
				}

				//  X
				// XY
				{
					auto& x1 = padded_grid[i + 1][j];
					auto& x2 = padded_grid[i][j + 1];
					auto& y = padded_grid[i + 1][j + 1];
					add_triple_constraints(x1, x2, y);
				}
			}
		}

		add_constraint(4, turn_sum, 4);

		// Forbid the two configurations
		//
		//  01     10
		//  10 and 01
		//
		// This is not completely correct, but the code
		// above can not handle them.
		for (int i = 0; i < m + 1; ++i) {
			for (int j = 0; j < n + 1; ++j) {
				add_constraint(padded_grid[i][i] + (1 - padded_grid[i + 1][j])
				                   + (1 - padded_grid[i][j + 1]) + padded_grid[i + 1][j + 1]
				               <= 3);
				add_constraint((1 - padded_grid[i][i]) + padded_grid[i + 1][j]
				                   + padded_grid[i][j + 1] + (1 - padded_grid[i + 1][j + 1])
				               <= 3);
			}
		}
	}

	// The second way is based on adding flow to one square and
	// requiring it to dissipate through the grid.
	void add_connected_constraint_flow(const vector<vector<BooleanVariable>>& grid) {
		ConnectedComponentFlow<> connected(this);
		int m = grid.size();
		int n = grid[0].size();
		for (int i : range(m)) {
			for (int j : range(n)) {
				connected.add_node(grid.at(i).at(j));
			}
		}

		for (int i : range(m)) {
			for (int j : range(n)) {
				if (i > 0) {
					connected.add_edge(grid[i - 1][j], grid[i][j]);
				}
				if (j > 0) {
					connected.add_edge(grid[i][j - 1], grid[i][j]);
				}
			}
		}

		vector<BooleanVariable> border_nodes;
		for (int i : range(m)) {
			border_nodes.push_back(grid[i][0]);
			border_nodes.push_back(grid[i][n - 1]);
		}
		for (int j : range(1, n - 1)) {
			border_nodes.push_back(grid[0][j]);
			border_nodes.push_back(grid[m - 1][j]);
		}

		connected.set_first_node(border_nodes);
		connected.finish();
	}

	void add_connected_constraint_noflow(const vector<vector<BooleanVariable>>& grid) {
		int m = grid.size();
		int n = grid[0].size();

		ConnectedComponentNoFlow<> connected(this, m * n / 2);

		for (int i : range(m)) {
			for (int j : range(n)) {
				connected.add_node(grid.at(i).at(j));
			}
		}

		for (int i : range(m)) {
			for (int j : range(n)) {
				if (i > 0) {
					connected.add_edge(grid[i - 1][j], grid[i][j]);
				}
				if (j > 0) {
					connected.add_edge(grid[i][j - 1], grid[i][j]);
				}
			}
		}

		connected.finish();
	}
};

int main_program() {
	using namespace std;

	int n = 6;

	MyIP ip;

	// Parking spots. We want to maximize the number of these.
	auto P = ip.add_boolean_grid(n, n, -1.0);
	// Transporation to parking spots.
	auto T = ip.add_boolean_grid(n, n);

	int starti = 0;
	int startj = 1;

	// Entrance to the parking lot.
	ip.add_constraint(T[starti][startj]);

	// Border
	// Sum border = 0;
	// for (int i = 1; i < n-1; ++i) {
	// 	border += T[0][i];
	// 	border += T[n-1][i];
	// 	border += T[i][0];
	// 	border += T[i][n-1];
	// }
	// // There has to be an entrance.
	// ip.add_constraint(border >= 1);

	// T[0][0] == 1
	//. P P P P P
	//. . . . . .
	//. P P P P P
	//. P P P P P
	//. . . . . .
	// P P P P P P
	// 21

	// T[0][1] == 1
	// P . P P P P
	// P . . . . .
	// P . P P P P
	// P . P P P P
	// P . . . . .
	//  P P P P P
	// 22

	// T[0][2] == 1
	// P . . . . P
	// P . P P . P
	// P . P P . P
	// P . P P . P
	// P . P P . P
	// P . P P . P
	// 22

	// T[0][1] == 1
	// P . P P P P P P . P
	// P . P . . . . . . P
	// P . P P P P P P . P
	// P . P   P P P P . P
	// P . P P . . . . . P
	// P . P P . P . P . P
	// P . P P . P P P . P
	// P . P P . P P P P P
	// . . . . . . . . . .
	// P P P P P P P P P P
	// 61

	// We can not have both parking and transportation.
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			// Not allowing both to be zero seems to give a slightly
			// stronger formulation.
			ip.add_constraint(P[i][j] + T[i][j] == 1);
		}
	}

	// Every parking lot needs to be adjacent to a transport square.
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Sum neighborhood;
			if (i > 0) {
				neighborhood += T[i - 1][j];
			}
			if (i < n - 1) {
				neighborhood += T[i + 1][j];
			}
			if (j > 0) {
				neighborhood += T[i][j - 1];
			}
			if (j < n - 1) {
				neighborhood += T[i][j + 1];
			}

			ip.add_constraint(P[i][j] <= neighborhood);
			ip.add_constraint(T[i][j] <= neighborhood + P[i][j]);
		}
	}

	// The hardest part: the transport squares need to be
	// connected to each other.
	// ip.add_connected_constraint_perimeter(T);  // Problem-specific.
	// ip.add_connected_constraint_flow(T);  // Best with an IP solver.
	ip.add_connected_constraint_noflow(T);  // Best with the SAT solver.

	auto print_solution = [&]() {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				auto p = ip.get_solution(P[i][j]);
				auto t = ip.get_solution(T[i][j]);
				if (p && t) {
					// Infeasible.
					cout << '#';
				} else if (!p && !t) {
					cout << ' ';
				} else if (p) {
					cout << 'P';
				} else {
					cout << '.';
				}
				cout << ' ';
			}
			cout << endl;
		}
		cout << endl;
	};

	IpToSatSolver solver(bind(minisat_solver, false));
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		cout << "\nSolution:\n";
		print_solution();

		break;
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
