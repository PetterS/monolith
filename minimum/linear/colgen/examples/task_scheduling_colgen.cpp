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

#include <minimum/algorithms/dag.h>
#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/proto.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

DEFINE_int32(num_solutions, 1, "Number of solutions to compute. Default: 1.");

DEFINE_bool(sat,
            false,
            "Use a SAT solver to try to obtain a feasible solution as soon as possible. The SAT "
            "problem is not easy to solve. Default: false.");

DEFINE_int32(verify_solution_value,
             -1,
             "If non-negative, verifies that some solution is better than this value. "
             "Default: -1.");

bool ok_together(const proto::TaskSchedulingProblem::Task& task1,
                 const proto::TaskSchedulingProblem::Task& task2) {
	return task2.start() > task1.end() || task1.start() > task2.end();
}
bool overlaps(const proto::TaskSchedulingProblem::Task& task1,
              const proto::TaskSchedulingProblem::Task& task2) {
	return !ok_together(task1, task2);
}

class TaskShedulingColgenProblem : public colgen::SetPartitioningProblem {
   public:
	TaskShedulingColgenProblem(const proto::TaskSchedulingProblem& problem_);

	virtual void generate(const std::vector<double>& dual_variables) override;

	virtual double integral_solution_value() override;

	virtual std::string member_name(int member) const override { return to_string('M', member); }

	std::vector<std::vector<int>> solution();

   private:
	const proto::TaskSchedulingProblem& problem;
};

TaskShedulingColgenProblem::TaskShedulingColgenProblem(const proto::TaskSchedulingProblem& problem_)
    : colgen::SetPartitioningProblem(problem_.machine_size(), problem_.task_size()),
      problem(problem_) {
	// Covering constraints with slack variables.
	for (auto c : range(problem.task_size())) {
		initialize_constraint(c, 1, 1, 0, 10'000);
	}

	for (auto m : range(problem.machine_size())) {
		colgen::Column column(0, 0, 1);
		column.add_coefficient(m, 1);
		add_column(move(column));
	}
}

void TaskShedulingColgenProblem::generate(const std::vector<double>& dual_variables) {
	vector<colgen::Column> columns_to_add(problem.machine_size());

#pragma omp parallel for
	for (int m = 0; m < problem.machine_size(); ++m) {
		int num_nodes = 0;
		int s = num_nodes++;
		vector<int> task_to_node(problem.task_size(), -1);
		for (int task : problem.machine(m).qualification()) {
			bool ok = true;
			for (int task2 : problem.machine(m).qualification()) {
				if (fixes_for_member(m).at(task2) == 1
				    && overlaps(problem.task(task), problem.task(task2))) {
					ok = false;
					break;
				}
			}
			if (fixes_for_member(m).at(task) == 0) {
				ok = false;
			}
			if (!ok) {
				continue;
			}

			task_to_node[task] = num_nodes++;
		}
		int t = num_nodes++;

		vector<int> node_to_task(num_nodes, -1);
		for (int task : problem.machine(m).qualification()) {
			int node = task_to_node[task];
			if (node >= 0) {
				node_to_task[task_to_node[task]] = task;
			}
		}
		minimum::algorithms::SortedDAG<double, 0> dag(num_nodes);

		for (int task : problem.machine(m).qualification()) {
			int node = task_to_node[task];
			if (node >= 0) {
				dag.add_edge(s, task_to_node[task]);
				dag.add_edge(task_to_node[task], t);
			}
		}
		for (int task1 : problem.machine(m).qualification()) {
			for (int task2 : problem.machine(m).qualification()) {
				int node1 = task_to_node[task1];
				int node2 = task_to_node[task2];
				if (node1 >= 0 && node2 >= 0
				    && problem.task(task1).start() < problem.task(task2).start()
				    && ok_together(problem.task(task1), problem.task(task2))) {
					dag.add_edge(task_to_node[task1], task_to_node[task2]);
				}
			}
		}
		for (int task : range(problem.task_size())) {
			int i = task_to_node[task];
			if (i >= 0) {
				dag.change_node_cost(i, -dual_variables[problem.machine_size() + task]);
			}
		}

		vector<minimum::algorithms::SolutionEntry<double>> solution;
		minimum::algorithms::shortest_path(dag, &solution);

		int i = dag.size() - 1;
		vector<int> tasks_to_add;
		tasks_to_add.reserve(solution.size());
		while (i >= 0) {
			int task = node_to_task[i];
			if (task >= 0) {
				tasks_to_add.push_back(task);
			}
			i = solution[i].prev;
		}

		// Make sure all tasks that are fixed are added. We could construct the
		// graph so that this is required, but that is a little more complicated.
		for (int task : problem.machine(m).qualification()) {
			if (fixes_for_member(m).at(task) == 1) {
				auto itr = find(tasks_to_add.begin(), tasks_to_add.end(), task);
				if (itr == tasks_to_add.end()) {
					tasks_to_add.push_back(task);
				}
			}
		}

		colgen::Column column(1, 0, 1);
		column.add_coefficient(m, 1);
		for (int task : tasks_to_add) {
			column.add_coefficient(problem.machine_size() + task, 1);
		}
		columns_to_add[m] = move(column);
	}

	for (int m = 0; m < problem.machine_size(); ++m) {
		add_column(move(columns_to_add[m]));
	}
}

double TaskShedulingColgenProblem::integral_solution_value() {
	if (!FLAGS_sat) {
		return SetPartitioningProblem::integral_solution_value();
	}

	IP ip;
	vector<Sum> column_to_sum(active_columns().size(), Sum(0));
	vector<Sum> machine_sums(problem.machine_size(), Sum(0));
	vector<Sum> task_sums(problem.task_size(), Sum(0));
	for (auto i : range(active_columns().size())) {
		auto col_index = active_columns()[i];
		auto m = column_member(col_index);
		if (m < 0) {
			continue;
		}
		auto x = ip.add_boolean();
		machine_sums[m] += x;
		column_to_sum[i] += x;

		auto& column = pool.at(col_index);
		if (column.is_fixed()) {
			ip.add_constraint(x == column.solution_value);
		}

		for (auto& entry : column) {
			if (entry.row >= problem.machine_size()) {
				int task = entry.row - problem.machine_size();
				task_sums[task] += x;
			}
		}
	}
	try {
		for (auto& sum : machine_sums) {
			ip.add_constraint(sum == 1);
		}
		for (auto& sum : task_sums) {
			ip.add_constraint(sum >= 1);
		}
	} catch (std::runtime_error& err) {
		return SetPartitioningProblem::integral_solution_value();
	}

	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	if (solver.solutions(&ip)->get()) {
		cerr << "-- SAT found a feasible solution.\n";
		double cost = 0;
		for (auto i : range(active_columns().size())) {
			if (column_to_sum[i].value() > 0.5) {
				pool.at(active_columns()[i]).fix(1);
				cost += pool.at(active_columns()[i]).cost();
			}
		}
		return cost;
	} else {
		return SetPartitioningProblem::integral_solution_value();
	}
}

vector<vector<int>> TaskShedulingColgenProblem::solution() {
	compute_fractional_solution(active_columns());
	vector<vector<int>> solution(problem.machine_size());

	for (auto m : range(problem.machine_size())) {
		for (int t : range(problem.task_size())) {
			if (fractional_solution(m, t) > 0.5) {
				solution[m].push_back(t);
			}
		}
	}

	return solution;
}

int main_program(int num_args, char* args[]) {
	using namespace std;
	check(num_args >= 2, "Need a file name.");
	ifstream fin(args[1]);
	auto problem = read_ptask_instance(fin);
	TaskShedulingColgenProblem colgen_problem(problem);

	vector<int> solution_objectives;
	for (int iter = 1; iter <= FLAGS_num_solutions; ++iter) {
		colgen_problem.unfix_all();
		colgen_problem.solve();

		auto solution = colgen_problem.solution();
		int machines_used = 0;
		set<int> tasks;
		for (auto m : range(problem.machine_size())) {
			cout << "Machine " << m << ": " << to_string(solution.at(m)) << "\n";
			if (!solution.at(m).empty()) {
				machines_used++;
			}
			for (int t : solution.at(m)) {
				tasks.insert(t);
			}
		}
		if (tasks.size() == problem.task_size()) {
			cout << "-- Machines used: " << machines_used << ".\n";
		} else {
			cout << "-- Infeasible.\n";
		}
		solution_objectives.push_back(machines_used);
	}

	cout << "-- All solutions: " << to_string(solution_objectives) << ".\n";

	if (FLAGS_verify_solution_value >= 0) {
		auto best = min_element(solution_objectives.begin(), solution_objectives.end());
		check(*best <= FLAGS_verify_solution_value, "Need a better solution.");
	}

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
