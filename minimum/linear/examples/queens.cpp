// Solves the n-queens problem.
//
// Try e.g. --n=8 --solver=constraint
#include <iostream>

#include <absl/container/flat_hash_map.h>
#include <gflags/gflags.h>

#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/range.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver_from_command_line.h>
using minimum::core::range;
using minimum::linear::Sum;

DEFINE_int32(n, 8, "Board size. Default: 8.");

int main_program(int num_args, char* args[]) {
	const int n = FLAGS_n;
	minimum::linear::IP ip;
	auto x = ip.add_boolean_grid(n, n);

	// Rows and columns.
	for (auto i : range(n)) {
		Sum row = 0;
		Sum col = 0;
		for (auto j : range(n)) {
			row += x[i][j];
			col += x[j][i];
		}
		ip.add_constraint(row == 1);
		ip.add_constraint(col == 1);
	}

	// Diagonals.
	absl::flat_hash_map<int, Sum> diagonals1;
	absl::flat_hash_map<int, Sum> diagonals2;
	for (auto i : range(n)) {
		for (auto j : range(n)) {
			diagonals1[i - j] += x[i][j];
			diagonals2[i + j] += x[i][j];
		}
	}
	for (const auto& entry : diagonals1) {
		ip.add_constraint(entry.second <= 1);
	}
	for (const auto& entry : diagonals2) {
		ip.add_constraint(entry.second <= 1);
	}

	auto solver = minimum::linear::create_solver_from_command_line();
	solver->set_silent(true);

	minimum::core::Timer timer("Solving");
	auto solutions = solver->solutions(&ip);
	int num_solutions = 0;
	while (solutions->get()) {
		num_solutions++;
		timer.check_for_interruption();
	}
	timer.OK();

	std::cout << num_solutions << " solutions found.\n";
	for (auto i : range(n)) {
		for (auto j : range(n)) {
			std::cout << (x[i][j].bool_value() ? "X " : ". ");
		}
		std::cout << "\n";
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
