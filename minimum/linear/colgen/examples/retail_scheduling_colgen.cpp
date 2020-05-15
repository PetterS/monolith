#include <fstream>
#include <random>

#include <gflags/gflags.h>

#include <minimum/algorithms/dag.h>
#include <minimum/core/color.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/random.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/retail_scheduling.h>
#include <minimum/linear/solver.h>

using namespace std;
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

DEFINE_int32(num_solutions, 1, "Number of solutions to compute. Default: 1.");

class GraphBuilder {
   public:
	GraphBuilder(const RetailProblem& problem_) : problem(problem_) {}

	int num_nodes() const { return 2 * problem.periods.size() + 2; }
	int source() const {
		minimum_core_assert(!is_not_worked(0) && !is_have_worked(0));
		return 0;
	}

	int sink() const {
		int node = num_nodes() - 1;
		minimum_core_assert(!is_not_worked(0) && !is_have_worked(0));
		return node;
	}

	int not_worked_node(int period) const {
		minimum_core_assert(0 <= period && period < problem.periods.size());
		int node = 1 + 2 * period;
		minimum_core_assert(get_period(node) == period);
		minimum_core_assert(is_not_worked(node) && !is_have_worked(node));
		return node;
	}

	int have_worked_node(int period) const {
		minimum_core_assert(0 <= period && period < problem.periods.size());
		int node = 1 + 2 * period + 1;
		minimum_core_assert(get_period(node) == period);
		minimum_core_assert(!is_not_worked(node) && is_have_worked(node));
		return node;
	}

	int get_period(int node) const {
		minimum_core_assert(1 <= node && node < num_nodes() - 1, "Node is ", node);
		return (node - 1) / 2;
	}

	bool is_not_worked(int node) const {
		return 1 <= node && node < num_nodes() - 1 && (node - 1) % 2 == 0;
	}

	bool is_have_worked(int node) const {
		return 1 <= node && node < num_nodes() - 1 && (node - 1) % 2 == 1;
	}

	// float may be OK if more efficient.
	minimum::algorithms::SortedDAG<double, 1, 1> build(const vector<double>& duals,
	                                                   int staff_index,
	                                                   const vector<vector<int>>& fixes,
	                                                   mt19937* rng) const {
		vector<int> best_task(problem.periods.size());
		vector<double> best_dual(problem.periods.size());
		for (auto p : range(problem.periods.size())) {
			best_dual[p] = numeric_limits<double>::min();
			for (auto t : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(p, t);
				int row = problem.staff.size() + c;
				double dual = duals.at(row);
				if (dual > best_dual[p]) {
					best_dual[p] = dual;
					best_task[p] = t;
				}
			}
		}

		vector<double> subset_sum(problem.periods.size());
		subset_sum[0] = best_dual[0];
		for (auto p : range(size_t{1}, problem.periods.size())) {
			subset_sum[p] = subset_sum[p - 1] + best_dual[p];
		}

		minimum::algorithms::SortedDAG<double, 1, 1> dag(num_nodes());

		dag.add_edge(source(), not_worked_node(0));
		dag.add_edge(not_worked_node(problem.periods.size() - 1), sink());

		for (int day = 0; day < problem.num_days; ++day) {
			int day_end = problem.day_start(day + 1);
			for (int p = problem.day_start(day); p < day_end; ++p) {
				// Edges from one "not worked" to the next "not worked" span days.
				if (p < problem.periods.size() - 1) {
					dag.add_edge(not_worked_node(p), not_worked_node(p + 1));
				}

				for (int shift_length = 6 * 4; shift_length <= 10 * 4; ++shift_length) {
					if (p + shift_length >= day_end) {
						break;
					}

					int last_node_in_shift = p + shift_length - 1;
					double shift_cost = subset_sum[last_node_in_shift];
					if (p > 0) {
						shift_cost -= subset_sum[p - 1];
					}
					auto& edge = dag.add_edge(
					    not_worked_node(p), have_worked_node(last_node_in_shift), shift_cost);
					edge.weights[0] = shift_length;
				}

				int next_day_in_14_hours = p + 14 * 4;
				if (next_day_in_14_hours < day_end) {
					next_day_in_14_hours = day_end;
				}
				int target_node;
				if (next_day_in_14_hours < problem.periods.size()) {
					target_node = not_worked_node(next_day_in_14_hours);
				} else {
					target_node = sink();
				}
				dag.add_edge(have_worked_node(p), target_node);
			}
		}

		return dag;
	}

   private:
	const RetailProblem& problem;
};

bool create_roster_graph(const RetailProblem& problem,
                         const vector<double>& duals,
                         int staff_index,
                         const vector<vector<int>>& fixes,
                         vector<vector<int>>* solution,
                         mt19937* rng) {
	GraphBuilder builder(problem);
	auto dag = builder.build(duals, staff_index, fixes, rng);

	vector<int> graph_solution;
	resource_constrained_shortest_path(dag,
	                                   problem.staff.at(staff_index).min_minutes / 15,
	                                   problem.staff.at(staff_index).max_minutes / 15,
	                                   &graph_solution);

	minimum_core_assert(problem.num_tasks == 1, "Need code for this.");
	for (int i = 0; i + 1 < graph_solution.size(); ++i) {
		int node1 = graph_solution[i];
		int node2 = graph_solution[i + 1];
		if (builder.is_not_worked(node1) && builder.is_have_worked(node2)) {
			int p1 = builder.get_period(node1);
			int p2 = builder.get_period(node2);
			for (int p = p1; p <= p2; ++p) {
				(*solution)[p][0] = 1;
			}
		}
	}
	return true;
}

bool create_roster_heuristic(const RetailProblem& problem,
                             const vector<double>& duals,
                             int staff_index,
                             const vector<vector<int>>& fixes,
                             vector<vector<int>>* solution,
                             mt19937* rng) {
	uniform_int_distribution random_day(0, problem.num_days - 1);

	for (auto p : range(problem.periods.size())) {
		for (auto t : range(problem.num_tasks)) {
			(*solution)[p][t] = 0;
		}
	}

	vector<int> best_task(problem.periods.size());
	vector<double> best_dual(problem.periods.size());
	vector<int> local_solution(problem.periods.size(), 0);

	for (auto p : range(problem.periods.size())) {
		best_dual[p] = numeric_limits<double>::min();
		for (auto t : range(problem.num_tasks)) {
			int c = problem.cover_constraint_index(p, t);
			int row = problem.staff.size() + c;
			double dual = duals.at(row);
			if (dual > best_dual[p]) {
				best_dual[p] = dual;
				best_task[p] = t;
			}
		}
	}

	vector<double> subset_sum(problem.periods.size());
	subset_sum[0] = best_dual[0];
	for (auto p : range(size_t{1}, problem.periods.size())) {
		subset_sum[p] = subset_sum[p - 1] + best_dual[p];
	}

	vector<int> shift_start(problem.num_days, -1);
	vector<int> shift_end(problem.num_days, -1);

	for (int day = 0; day < problem.num_days; ++day) {
		double best_value = numeric_limits<double>::min();
		int best_start = -1;
		int best_end = -1;

		int initial_p = problem.day_start(day);
		if (day > 0 && shift_end[day - 1] >= 0) {
			initial_p = max(initial_p, shift_end[day - 1] + 14 * 4 + 1);
		}

		for (int p_start = initial_p;
		     p_start < problem.day_start(day + 1) && problem.periods[p_start].day == day;
		     ++p_start) {
			if (!problem.periods[p_start].can_start) {
				continue;
			}

			for (int p_end = problem.day_start(day) + problem.min_shift_length - 1;
			     p_end < problem.day_start(day + 1)
			     && (p_end - p_start + 1) <= problem.max_shift_length;
			     ++p_end) {
				double value = subset_sum[p_end];
				if (p_start > 0) {
					value -= subset_sum[p_start - 1];
				}

				if (value > best_value) {
					best_value = value;
					best_start = p_start;
					best_end = p_end;
				}
			}
		}

		shift_start[day] = best_start;
		shift_end[day] = best_end;
	}

	for (int iteration = 1; iteration <= 100; iteration++) {
		bool feasible = true;

		int working_minutes = 0;
		for (int day : range(problem.num_days)) {
			if (shift_start[day] >= 0) {
				working_minutes += (shift_end[day] - shift_start[day] + 1) * 15;
			}
		}

		if (working_minutes < problem.staff[staff_index].min_minutes) {
			feasible = false;

			// TODO: Go through the days in random order.
			for (int day : range(problem.num_days)) {
				// TODO: Will violate shift length constraints.
				if (shift_start[day] > problem.day_start(day)) {
					shift_start[day]--;
				}
				if (shift_end[day] < problem.day_start(day + 1)) {
					shift_end[day]++;
				}

				if (working_minutes >= problem.staff[staff_index].min_minutes) {
					break;
				}
			}
		} else if (working_minutes > problem.staff[staff_index].max_minutes) {
			feasible = false;
			int delta = working_minutes - problem.staff[staff_index].max_minutes;

			int day = random_day(*rng);
			shift_start[day] = -1;
			shift_end[day] = -2;
		}

		if (feasible) {
			for (int day : range(problem.num_days)) {
				if (shift_start[day] >= 0) {
					for (int p = shift_start[day]; p < shift_end[day]; ++p) {
						(*solution)[p][best_task[p]] = 1;
					}
				}
			}
			return true;
		}
	}

	return false;
}

bool create_roster_dummy(const RetailProblem& problem,
                         const vector<double>& duals,
                         int staff_index,
                         const vector<vector<int>>& fixes,
                         vector<vector<int>>* solution,
                         mt19937* rng) {
	uniform_real_distribution<double> rand(0, 1);

	// TODO: Implement this method with something reasonable.

	for (auto d : range(problem.periods.size())) {
		for (auto s : range(problem.num_tasks)) {
			if (fixes[d][s] >= 0) {
				(*solution)[d][s] = fixes[d][s];
				continue;
			}

			int c = problem.cover_constraint_index(d, s);
			int row = problem.staff.size() + c;
			if (duals.at(row) > 0) {
				(*solution)[d][s] = 1;
			} else {
				(*solution)[d][s] = 0;
			}
		}
	}

	return true;
}

bool create_roster_ip(const RetailProblem& problem,
                      const vector<double>& duals,
                      int staff_index,
                      const vector<vector<int>>& fixes,
                      vector<vector<int>>* solution,
                      mt19937* rng) {
	uniform_real_distribution<double> rand(0, 1);

	IP ip;
	auto working = ip.add_boolean_grid(problem.periods.size(), problem.num_tasks);

	for (int p = 0; p < problem.periods.size(); ++p) {
		Sum task_sum = 0;
		for (int t = 0; t < problem.num_tasks; ++t) {
			task_sum += working[p][t];
		}
		ip.add_constraint(task_sum <= 1);
	}

	vector<Sum> working_any_shift(problem.periods.size(), 0);
	for (int p = 0; p < problem.periods.size(); ++p) {
		for (int t = 0; t < problem.num_tasks; ++t) {
			working_any_shift[p] += working[p][t];
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

	for (int p = 0; p < problem.periods.size(); ++p) {
		for (int t = 0; t < problem.num_tasks; ++t) {
			if (fixes[p][t] >= 0) {
				ip.add_constraint(working[p][t] == fixes[p][t]);
			} else {
				int c = problem.cover_constraint_index(p, t);
				int row = problem.staff.size() + c;
				ip.add_objective(-duals.at(row) * working[p][t]);
			}
		}
	}

	IPSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	check(solutions->get(), "Could not solve problem.");

	for (int p = 0; p < problem.periods.size(); ++p) {
		for (int t = 0; t < problem.num_tasks; ++t) {
			(*solution)[p][t] = working[p][t].bool_value() ? 1 : 0;
		}
	}

	return true;
}

class ShiftShedulingColgenProblem : public SetPartitioningProblem {
   public:
	ShiftShedulingColgenProblem(RetailProblem problem_)
	    : SetPartitioningProblem(problem_.staff.size(), problem_.num_cover_constraints()),
	      problem(move(problem_)) {
		for (int p : range(problem.staff.size())) {
			rng.emplace_back(repeatably_seeded_engine<std::mt19937>(p));
		}

		Timer timer("Setting up colgen problem");
		solution = make_grid<int>(problem.staff.size(), problem.periods.size(), problem.num_tasks);

		// Covering constraints with slack variables.
		for (auto d : range(problem.periods.size())) {
			auto& period = problem.periods[d];
			for (auto s : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(d, s);
				// Initialize this constraint.
				initialize_constraint(
				    c, period.min_cover.at(s), period.max_cover.at(s), 1000.0, 1.0);
			}
		}
		timer.OK();

		timer.next("Creating initial rosters");

		// Initial duals just try to cover as much as possible, with some
		// randomness to encourage different rosters.
		auto initial_duals = make_grid<double>(problem.staff.size(), number_of_rows());
		uniform_real_distribution<double> dual(-100, 100);
		for (int p = 0; p < problem.staff.size(); ++p) {
			for (int r = 0; r < number_of_rows(); ++r) {
				initial_duals[p][r] = dual(rng[p]);
			}
		}

		auto no_fixes =
		    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });

		for (int p = 0; p < problem.staff.size(); ++p) {
			minimum_core_assert(
			    create_roster_graph(problem, initial_duals[p], p, no_fixes, &solution[p], &rng[p]));
		}

		for (int p = 0; p < problem.staff.size(); ++p) {
			add_column(create_column(p, solution[p]));
		}

		timer.OK();
	}

	~ShiftShedulingColgenProblem() {}

	Column create_column(int staff_index, const vector<vector<int>>& solution_for_staff) {
		double cost = 0;

		// TODO: Set the cost correctly.

		Column column(cost, 0, 1);
		column.add_coefficient(staff_index, 1);

		int start_row = problem.staff.size();
		for (auto d : range(problem.periods.size())) {
			for (auto s : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(d, s);
				int sol = solution_for_staff[d][s];
				minimum_core_assert(sol == 0 || sol == 1);
				if (sol == 1) {
					column.add_coefficient(start_row + c, 1.0);
				}
			}
		}
		return column;
	}

	virtual void generate(const std::vector<double>& dual_variables) override {
		// vector of ints instead of bool because we may be writing to
		// it concurrently.
		vector<int> success(problem.staff.size(), 0);

		for (int p = 0; p < problem.staff.size(); ++p) {
			if (generate_for_staff(p, dual_variables)) {
				success[p] = 1;
			}
		}

		for (int p = 0; p < problem.staff.size(); ++p) {
			if (success[p] == 1) {
				add_column(create_column(p, solution[p]));
			}
		}
	}

	virtual std::string member_name(int member) const override {
		return problem.staff.at(member).id;
	}

	const std::vector<std::vector<std::vector<int>>>& get_solution() const { return solution; }

   private:
	bool generate_for_staff(int staff_index, const std::vector<double>& dual_variables) {
		if (member_fully_fixed(staff_index)) {
			return false;
		}

		auto fixes = make_grid<int>(problem.periods.size(), problem.num_tasks);
		auto& fixes_for_staff = fixes_for_member(staff_index);
		for (auto c : range(problem.num_cover_constraints())) {
			int d = problem.cover_constraint_to_period(c);
			int s = problem.cover_constraint_to_task(c);
			fixes[d][s] = fixes_for_staff[c];
		}

		return create_roster_graph(
		    problem, dual_variables, staff_index, fixes, &solution[staff_index], &rng[staff_index]);
	}

	const minimum::linear::RetailProblem problem;
	std::vector<std::vector<std::vector<int>>> solution;
	vector<mt19937> rng;
};

int main_program(int num_args, char* args[]) {
	if (num_args <= 1) {
		cerr << "Need a file name.\n";
		return 1;
	}
	string filename_base = args[1];
	ifstream input(filename_base + ".txt");
	const RetailProblem problem(input);
	problem.print_info();

	ShiftShedulingColgenProblem colgen_problem(problem);

	auto start_time = wall_time();
	for (int i = 1; i <= FLAGS_num_solutions; ++i) {
		colgen_problem.unfix_all();
		// colgen_problem.solve();
		cerr << "Colgen done.\n";
		auto elapsed_time = wall_time() - start_time;
		cerr << "Elapsed time : " << elapsed_time << "s.\n";
	}

	problem.save_solution(
	    filename_base + ".ros", filename_base + ".xml", colgen_problem.get_solution());

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
