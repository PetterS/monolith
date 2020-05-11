#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <string_view>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/main.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/retail_scheduling.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/solver.h>

using namespace minimum::core;
using namespace minimum::linear;

int main_program(int num_args, char* args[]) {
	using namespace std;
	if (num_args <= 1) {
		cerr << "Need a file name.\n";
		return 1;
	}
	ifstream input(args[1]);
	const RetailProblem problem(input);
	problem.print_info();

	IP ip;
	auto working =
	    ip.add_boolean_cube(problem.staff.size(), problem.periods.size(), problem.num_tasks);
	for (int s = 0; s < problem.staff.size(); ++s) {
		for (int p = 0; p < problem.periods.size(); ++p) {
			Sum task_sum = 0;
			for (int t = 0; t < problem.num_tasks; ++t) {
				task_sum += working[s][p][t];
			}
			ip.add_constraint(task_sum <= 1);
		}
	}

	for (int s = 0; s < problem.staff.size(); ++s) {
		vector<Sum> working_any_shift(problem.periods.size(), 0);
		for (int p = 0; p < problem.periods.size(); ++p) {
			for (int t = 0; t < problem.num_tasks; ++t) {
				working_any_shift[p] += working[s][p][t];
			}
		}

		auto starting_shift = ip.add_boolean_vector(problem.periods.size());
		for (int p = 1; p < problem.periods.size(); ++p) {
			const auto& yesterday = working_any_shift[p - 1];
			const auto& today = working_any_shift[p];
			// today ∧ ¬yesterday ⇔ starting_shift[p]
			ip.add_constraint(today >= starting_shift[p]);
			ip.add_constraint(1 - yesterday >= starting_shift[p]);
			ip.add_constraint(starting_shift[p] + 1 >= today + (1 - yesterday));
		}

		vector<Sum> day_constraints(problem.num_days + 1, 0);
		for (int p = 1; p < problem.periods.size(); ++p) {
			auto& period = problem.periods.at(p);
			day_constraints[period.day] += starting_shift[p];
			if (!period.can_start) {
				ip.add_constraint(starting_shift[p] == 0);
			}
		}

		// Can only start a single shift each day.
		for (int day = 0; day < problem.num_days; ++day) {
			ip.add_constraint(day_constraints[day] <= 1);
		}
		ip.add_constraint(day_constraints.back() == 0);

		// Minimum shift duration is 6 hours.
		for (int p = 1; p < problem.periods.size(); ++p) {
			for (int p2 = p + 1; p2 < p + 6 * 4 && p2 < problem.periods.size(); ++p2) {
				ip.add_constraint(working_any_shift[p2] >= starting_shift[p]);
			}
		}
	}

	Sum objective;
	for (int p = 0; p < problem.periods.size(); ++p) {
		auto& period = problem.periods.at(p);
		for (int t = 0; t < problem.num_tasks; ++t) {
			Sum num_workers = 0;
			for (int s = 0; s < problem.staff.size(); ++s) {
				num_workers += working[s][p][t];
			}
			auto slack = ip.add_variable(IP::Real);
			ip.add_constraint(slack >= 0);
			ip.add_objective(slack);
			ip.add_constraint(period.min_cover[t] <= num_workers);
			ip.add_constraint(num_workers == period.max_cover[t] + slack);
		}
	}

	IPSolver solver;
	clog << "Starting search." << endl;
	auto solutions = solver.solutions(&ip);
	if (solutions->get()) {
		clog << "Objective " << objective.value() << " proven optimal.\n";
	} else {
		clog << "IP was not solved." << endl;
		return 2;
	}

	for (int p = 0; p < problem.periods.size(); ++p) {
		auto& period = problem.periods.at(p);
		cout << period.human_readable_time << ": ";
		for (int s = 0; s < problem.staff.size(); ++s) {
			bool printed = false;
			for (int t = 0; t < problem.num_tasks; ++t) {
				if (working[s][p][t].bool_value()) {
					minimum_core_assert(!printed);
					cout << t;
					printed = true;
				}
			}
			if (!printed) {
				cout << ".";
			}
		}
		cout << "\n";
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
