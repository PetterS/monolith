//
// A solver for Hashiwokakero puzzles.
// http://en.wikipedia.org/wiki/Hashiwokakero
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	const int _ = -1;
	const vector<vector<int>> given = {
	    {2, _, 2, _, 2, _, _, 2, _, 2, _, _, 2, _, _, _, _, 2, _, 2, _, 2, _, 2, _},
	    {_, 1, _, _, _, _, 2, _, _, _, 4, _, _, 5, _, 2, _, _, 1, _, 2, _, 2, _, 1},
	    {2, _, _, 5, _, 4, _, _, 3, _, _, _, _, _, 1, _, _, 4, _, 5, _, 1, _, 1, _},
	    {_, _, _, _, _, _, _, _, _, _, _, 1, _, 3, _, _, 1, _, _, _, _, _, _, _, _},
	    {2, _, _, 6, _, 6, _, _, 8, _, 5, _, 2, _, _, 3, _, 5, _, 7, _, _, 2, _, _},
	    {_, 1, _, _, _, _, _, _, _, _, _, 1, _, _, 2, _, _, _, _, _, 1, _, _, _, 3},
	    {2, _, _, _, _, 5, _, _, 6, _, 4, _, _, 2, _, _, _, 2, _, 5, _, 4, _, 2, _},
	    {_, 2, _, 2, _, _, _, _, _, _, _, _, _, _, _, 3, _, _, 3, _, _, _, 1, _, 2},
	    {_, _, _, _, _, _, _, _, _, _, 4, _, 2, _, 2, _, _, 1, _, _, _, 3, _, 1, _},
	    {2, _, 3, _, _, 6, _, _, 2, _, _, _, _, _, _, _, _, _, _, 3, _, _, _, _, _},
	    {_, _, _, _, 1, _, _, 2, _, _, 5, _, _, 1, _, 4, _, 3, _, _, _, _, 2, _, 4},
	    {_, _, 2, _, _, 1, _, _, _, _, _, _, 5, _, 4, _, _, _, _, 4, _, 3, _, _, _},
	    {2, _, _, _, 3, _, 1, _, _, _, _, _, _, _, _, 3, _, _, 5, _, 5, _, _, 2, _},
	    {_, _, _, _, _, 2, _, 5, _, _, 7, _, 5, _, 3, _, 1, _, _, 1, _, _, 1, _, 4},
	    {2, _, 5, _, 3, _, _, _, _, 1, _, 2, _, 1, _, _, _, _, 2, _, 4, _, _, 2, _},
	    {_, _, _, _, _, 1, _, _, _, _, _, _, _, _, _, _, 2, _, _, 2, _, 1, _, _, 3},
	    {2, _, 6, _, 6, _, _, 2, _, _, 2, _, 2, _, 5, _, _, _, _, _, 2, _, _, _, _},
	    {_, _, _, _, _, 1, _, _, _, 3, _, _, _, _, _, 1, _, _, 1, _, _, 4, _, 3, _},
	    {_, _, 4, _, 5, _, _, 2, _, _, _, 2, _, _, 6, _, 6, _, _, 3, _, _, _, _, 3},
	    {2, _, _, _, _, _, _, _, _, _, 2, _, _, 1, _, _, _, _, _, _, 1, _, _, 1, _},
	    {_, _, 3, _, _, 3, _, 5, _, 5, _, _, 4, _, 6, _, 7, _, _, 4, _, 6, _, _, 4},
	    {2, _, _, _, 3, _, 5, _, 2, _, 1, _, _, _, _, _, _, _, _, _, _, _, _, _, _},
	    {_, _, _, _, _, _, _, _, _, 1, _, _, _, _, _, _, 3, _, 2, _, _, 5, _, _, 5},
	    {2, _, 3, _, 3, _, 5, _, 4, _, 3, _, 3, _, 4, _, _, 2, _, 2, _, _, _, 1, _},
	    {_, 1, _, 2, _, 2, _, _, _, 2, _, 2, _, _, _, 2, _, _, _, _, 2, _, 2, _, 2}};

	const int m = given.size();
	const int n = given.back().size();

	struct Edge {
		int i1;
		int i2;
		int j1;
		int j2;

		bool has_endpoint(int i, int j) const {
			return (i == i1 && j == j1) || (i == i2 && j == j2);
		}

		bool contains_point(int i, int j) const { return i1 <= i && i <= i2 && j1 <= j && j <= j2; }

		bool vertical() const { return j1 == j2; }
	};

	IP ip;
	vector<Edge> edges;
	vector<BooleanVariable> single_edge;
	vector<BooleanVariable> double_edge;

	// Generate all possible edges.
	int num_givens = 0;
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] == _) {
				continue;
			}

			num_givens++;

			for (int i2 = i + 1; i2 < m; i2++) {
				if (given[i2][j] != _) {
					edges.push_back({i, i2, j, j});
					single_edge.emplace_back(ip.add_boolean());
					double_edge.emplace_back(ip.add_boolean());
					break;
				}
			}

			for (int j2 = j + 1; j2 < n; j2++) {
				if (given[i][j2] != _) {
					edges.push_back({i, i, j, j2});
					single_edge.emplace_back(ip.add_boolean());
					double_edge.emplace_back(ip.add_boolean());
					break;
				}
			}
		}
	}
	auto num_edges = edges.size();
	minimum_core_assert(num_edges == single_edge.size());
	minimum_core_assert(num_edges == double_edge.size());

	//
	// Use a single edge or a double edge.
	//
	for (auto e : range(num_edges)) {
		ip.add_constraint(single_edge[e] + double_edge[e] <= 1);
	}

	//
	// The edges used should match what is given.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] == _) {
				continue;
			}

			Sum incoming = 0;
			for (auto e : range(num_edges)) {
				if (edges[e].has_endpoint(i, j)) {
					incoming += single_edge[e] + 2 * double_edge[e];
				}
			}
			ip.add_constraint(incoming == given[i][j]);
		}
	}

	//
	// No crossings.
	//
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] != _) {
				continue;
			}

			Sum crossing_edges = 0;
			for (auto e : range(num_edges)) {
				if (edges[e].contains_point(i, j)) {
					crossing_edges += single_edge[e] + double_edge[e];
				}
			}
			ip.add_constraint(crossing_edges <= 1);
		}
	}

	//
	// Connectivity
	//
	ConnectedComponentFlow<pair<int, int>> connected(&ip);
	for (int i : range(m)) {
		for (int j : range(n)) {
			if (given[i][j] == _) {
				continue;
			}
			connected.add_node({i, j});
		}
	}
	for (auto e : range(num_edges)) {
		connected.add_edge({edges[e].i1, edges[e].j1},
		                   {edges[e].i2, edges[e].j2},
		                   single_edge[e] || double_edge[e]);
	}
	connected.finish();

	IPSolver solver;
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		cout << "\nSolution:\n";
		for (int i : range(m)) {
			for (int j : range(n)) {
				string value = "  ";
				if (given[i][j] != _) {
					value = " " + to_string(given[i][j]);
				} else {
					for (auto e : range(num_edges)) {
						if (edges[e].contains_point(i, j)) {
							if (single_edge[e].bool_value()) {
								value = edges[e].vertical() ? " |" : "--";
							} else if (double_edge[e].bool_value()) {
								value = edges[e].vertical() ? " \"" : "==";
							}
						}
					}
				}
				cout << value;
			}
			cout << "\n";
		}

		// break to save time.
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
