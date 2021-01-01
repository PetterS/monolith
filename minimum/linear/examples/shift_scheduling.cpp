// Petter Strandmark.
//
// Reads and solves problem files from “Employee Shift Scheduling Benchmark
// Data Sets”
// http://www.cs.nott.ac.uk/~tec/NRP/
//
// First argument to the program is the problem .txt file.
//
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

std::mt19937_64 rng;

DEFINE_bool(save_mps, false, "Save an MPS file with the IP.");
DEFINE_bool(first_order_solver, false, "Use a first-order solver for the relaxation.");

int main_program(int num_args, char* args[]) {
	using namespace std;

	check(num_args == 2, "Need a file name.");
	ifstream fin(args[1]);
	auto problem = read_nottingham_instance(fin);

	double start_time = minimum::core::wall_time();

	IP ip;

	// The assignment variables.
	//   x[p][d][s] = 1 ⇔ Person p is working day d, shift s.
	auto x = ip.add_boolean_cube(problem.worker_size(), problem.num_days(), problem.shift_size());
	auto objective = create_ip(ip, problem, x);
	IPSolver solver;

	if (FLAGS_save_mps) {
		Timer t("Saving MPS");
		solver.save_MPS(ip, args[1]);  // .mps will be appended.
		t.OK();
		return 0;
	}

	bool first_solution = true;
	auto print_solution = [&]() {
		minimum_core_assert(ip.is_feasible_and_integral());
		double elapsed_time = minimum::core::wall_time() - start_time;
		clog << "Integer solution with objective " << objective.value() << " in " << elapsed_time
		     << " seconds." << endl;

		if (!first_solution) {
			cout << "======" << endl;
		}
		first_solution = false;

		for (int p = 0; p < problem.worker_size(); ++p) {
			for (int d = 0; d < problem.num_days(); ++d) {
				for (int s = 0; s < problem.shift_size(); ++s) {
					if (x[p][d][s].bool_value()) {
						cout << problem.shift(s).id();
					}
				}
				cout << '\t';
			}
			cout << endl;
		}
	};

	if (FLAGS_first_order_solver) {
		ScsSolver scs;
		auto solutions = scs.solutions(&ip);
		if (solutions->get()) {
			double elapsed_time = minimum::core::wall_time() - start_time;
			clog << "Objective " << objective.value() << " proven optimal in " << elapsed_time
			     << " seconds." << endl;
		} else {
			clog << "IP was not solved." << endl;
			return 2;
		}
		return 0;
	}

	minimum_core_assert(solver.solve_relaxation(&ip));
	clog << endl << "Lower bound is " << objective.value() << endl << endl;

	clog << "Starting search." << endl;

	auto solutions = solver.solutions(&ip);
	if (solutions->get()) {
		double elapsed_time = minimum::core::wall_time() - start_time;
		clog << "Objective " << objective.value() << " proven optimal in " << elapsed_time
		     << " seconds." << endl;
	} else {
		clog << "IP was not solved." << endl;
		return 2;
	}

	print_solution();

	return 0;
}

int main(int num_args, char* args[]) {
	if (num_args == 1 || ends_with(to_string(args[1]), "-help")) {
		gflags::SetUsageMessage(
		    "Runs direct IP method for shift scheduling.\n   shift_scheduling [options] "
		    "<input_file>");
		gflags::ShowUsageWithFlagsRestrict("shift_scheduling", "colgen");
		return 0;
	}
	return main_runner(main_program, num_args, args);
}
