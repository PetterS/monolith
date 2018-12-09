// Solves the n-queens domination problem.
// https://ajc.maths.uq.edu.au/pdf/15/ocr-ajc-v15-p145.pdf
// http://yetanothermathprogrammingconsultant.blogspot.com/2018/12/more-queens-problems.html
//
//
// Try e.g. --n=8 --solver=minisat
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
		for (auto j : range(n)) {
			ip.add_objective(x[i][j]);

			Sum sum = 0;
			auto add = [&sum, &x, n](int i2, int j2) {
				if (0 <= i2 && i2 < n && 0 <= j2 && j2 < n) {
					sum += x[i2][j2];
				}
			};
			for (auto k : range(n)) {
				add(k, j);
				add(i, k);
				add(i - k, j - k);
				add(i + k, j + k);
				add(i - k, j + k);
				add(i + k, j - k);
			}
			ip.add_constraint(sum >= 1);
		}
	}

	auto solver = minimum::linear::create_solver_from_command_line();
	solver->set_silent(true);

	minimum::core::Timer timer("Solving once");
	auto solutions = solver->solutions(&ip);
	minimum::core::check(solutions->get());
	timer.OK();

	for (auto i : range(n)) {
		for (auto j : range(n)) {
			std::cout << (x[i][j].bool_value() ? "X " : ". ");
		}
		std::cout << "\n";
	}

	timer.next("Finding all solutions");
	int num_solutions = 1;
	while (solutions->get()) {
		num_solutions++;
		timer.check_for_interruption();
	}
	timer.OK();
	std::cout << num_solutions << " solutions found.\n";
	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
