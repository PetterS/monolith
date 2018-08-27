//
// A solver for Polyomino puzzles.
// https://stackoverflow.com/questions/47918792/2d-bin-packing-on-a-grid
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/range.h>
#include <minimum/core/text_table.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program(int num_args, char* args[]) {
	using namespace std;

	constexpr int n = 29;

	// clang-format off
	constexpr int X = 1;
	constexpr int _ = 0;
	const vector<vector<vector<int>>> pieces = {
		{
			{ X },
			{ X },
			{ X },
			{ X }
		},
		{
			{ _, X, X},
			{ X, X, _ },
			{ _, X, _ }
		},
		{
			{ X, X, _ },
			{ _, X, X },
			{ _, X, _ }
		},
		{
			{ X, _ },
			{ X, _ },
			{ X, _ },
			{ X, X }
		},
		{
			{ _, X },
			{ _, X },
			{ _, X },
			{ X, X }
		},
		{
			{ X, X },
			{ X, X },
			{ _, X }
		},
		{
			{ X, X },
			{ X, X },
			{ X, _ }
		},
		{
			{ _, X },
			{ _, X },
			{ X, X },
			{ X, _ }
		},
		{
			{ X, _ },
			{ X, _ },
			{ X, X },
			{ _, X }
		},
	};
	// clang-format on

	IP ip;
	auto x = ip.add_boolean_cube(pieces.size(), n, n);
	auto grid = make_grid<Sum>(n, n);
	for (auto p : range(pieces.size())) {
		for (auto i : range(n)) {
			for (auto j : range(n)) {
				if (i + pieces[p].size() > n) {
					continue;
				}
				if (j + pieces[p][0].size() > n) {
					continue;
				}

				for (auto di : range(pieces[p].size())) {
					for (auto dj : range(pieces[p][0].size())) {
						if (pieces[p][di][dj] == X) {
							grid[i + di][j + dj] += x[p][i][j];
						}
					}
				}
			}
		}
	}

	for (auto i : range(n)) {
		for (auto j : range(n)) {
			ip.add_constraint(grid[i][j] == 1);
		}
	}

	IPSolver ip_solver;
	IpToSatSolver sat_solver([]() {
		auto solver = minisat_solver(false);
		// We don't need to use simplification since this problem
		// is pretty much optimally encoded for SAT.
		// Simplification will eliminate nothing.
		solver->use_simplification(false);
		return solver;
	});
	auto solutions = sat_solver.solutions(&ip);
	check(solutions->get(), "Solution failed.");

	auto solution = make_grid<int>(n, n);
	for (auto p : range(pieces.size())) {
		for (auto i : range(n)) {
			for (auto j : range(n)) {
				if (i + pieces[p].size() > n) {
					continue;
				}
				if (j + pieces[p][0].size() > n) {
					continue;
				}

				for (auto di : range(pieces[p].size())) {
					for (auto dj : range(pieces[p][0].size())) {
						if (pieces[p][di][dj] == X) {
							solution[i + di][j + dj] += (p + 1) * x[p][i][j].value();
						}
					}
				}
			}
		}
	}

	cout << "\nSolution:\n";
	print_table(cout, solution);

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
