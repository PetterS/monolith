// Petter Strandmark 2014.
//
// Reads and solves problem files from “Employee Shift Scheduling Benchmark
// Data Sets”
// http://www.cs.nott.ac.uk/~tec/NRP/
//
// First argument to the program is the problem .dat file.
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
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/proto.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

bool ok_together(const proto::TaskSchedulingProblem::Task& task1,
                 const proto::TaskSchedulingProblem::Task& task2) {
	return task2.start() > task1.end() || task1.start() > task2.end();
}
bool overlaps(const proto::TaskSchedulingProblem::Task& task1,
              const proto::TaskSchedulingProblem::Task& task2) {
	return !ok_together(task1, task2);
}

int main_program(int num_args, char* args[]) {
	using namespace std;
	check(num_args >= 2, "Need a file name.");
	ifstream fin(args[1]);
	auto problem = read_ptask_instance(fin);

	IP ip;
	auto x = ip.add_boolean_grid(problem.machine_size(), problem.task_size());

	// Machine qualifications.
	for (int m : range(problem.machine_size())) {
		for (int t : range(problem.task_size())) {
			bool found = false;
			for (int t2 : problem.machine(m).qualification()) {
				if (t == t2) {
					found = true;
					break;
				}
			}
			if (!found) {
				ip.add_constraint(x[m][t] == 0);
			}
		}
	}

	// Tasks can not overlap on the same machine.
	for (int m : range(problem.machine_size())) {
		auto used = ip.add_boolean();
		for (int t1 : range(problem.task_size())) {
			for (int t2 : range(problem.task_size())) {
				if (t1 != t2 && overlaps(problem.task(t1), problem.task(t2))) {
					ip.add_constraint(x[m][t1] + x[m][t2] <= used);
				}
			}
		}
		ip.add_objective(used);
	}

	// All tasks must be performed.
	for (int t : range(problem.task_size())) {
		Sum task_sum = 0;
		for (int m : range(problem.machine_size())) {
			task_sum += x[m][t];
		}
		ip.add_constraint(task_sum == 1);
	}

	IPSolver solver;
	check(solver.solutions(&ip)->get(), "Could not solve ip.");

	int machines_used = 0;
	for (int m : range(problem.machine_size())) {
		cout << "Machine " << m << ": ";
		bool used = false;
		for (int t : range(problem.task_size())) {
			if (x[m][t].bool_value()) {
				cout << t << " ";
				used = true;
			}
		}
		if (used) {
			machines_used++;
		}
		cout << "\n";
	}
	cout << "-- " << machines_used << " machines used.\n";

	return 0;
}

int main(int num_args, char* args[]) {
	gflags::ParseCommandLineFlags(&num_args, &args, true);
	try {
		return main_program(num_args, args);
	} catch (std::exception& err) {
		cerr << "Error: " << err.what() << endl;
		return 1;
	}
}
