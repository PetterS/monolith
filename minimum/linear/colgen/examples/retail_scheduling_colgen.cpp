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
		uniform_real_distribution<double> eps(-1, 1);

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
				auto& period = problem.periods.at(p);

				// Edges from one "not worked" to the next "not worked" span days.
				if (p < problem.periods.size() - 1) {
					dag.add_edge(not_worked_node(p), not_worked_node(p + 1));
				}

				// Shifts can only start between the following times: 06:00-10:00, 14:00-18:00 and
				// 20:00-00:00
				int start_on_day_hour4 = problem.time_to_hour4(period.human_readable_time);
				if ((6 * 4 <= start_on_day_hour4 && start_on_day_hour4 <= 10 * 4)
				    || (14 * 4 <= start_on_day_hour4 && start_on_day_hour4 <= 18 * 4)
				    || (20 * 4 <= start_on_day_hour4 && start_on_day_hour4 <= 24 * 4)) {
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
						    not_worked_node(p), have_worked_node(last_node_in_shift), -shift_cost);
						edge.weights[0] = shift_length;
					}

					int next_day_in_14_hours = p + 14 * 4 + 1;
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
	for (auto p : range(problem.periods.size())) {
		for (auto t : range(problem.num_tasks)) {
			(*solution)[p][t] = 0;
		}
	}

	GraphBuilder builder(problem);
	auto dag = builder.build(duals, staff_index, fixes, rng);

	vector<int> graph_solution;
	bool feasible = false;
	for (int attempts = 0; attempts <= 10; ++attempts) {
		resource_constrained_shortest_path(dag,
		                                   problem.staff.at(staff_index).min_minutes / 15,
		                                   problem.staff.at(staff_index).max_minutes / 15,
		                                   &graph_solution);
		feasible = true;

		//
		// Check consecutive working days.
		//
		vector<int> working_on_day(problem.num_days + 1, 0);
		for (int i = 0; i + 1 < graph_solution.size(); ++i) {
			int node1 = graph_solution[i];
			int node2 = graph_solution[i + 1];
			if (builder.is_not_worked(node1) && builder.is_have_worked(node2)) {
				int p1 = builder.get_period(node1);
				int day = problem.periods.at(p1).day;
				working_on_day[day] = 1;
			}
		}
		int consecutive = 0;
		vector<int> consecutive_working_days;
		consecutive_working_days.reserve(problem.num_days);
		for (int day : range(problem.num_days + 1)) {
			if (working_on_day[day] == 1) {
				consecutive_working_days.push_back(day);
			} else {
				if (consecutive_working_days.size() > 5) {
					feasible = false;
					uniform_int_distribution<int> random_index(0,
					                                           consecutive_working_days.size() - 1);
					int day = consecutive_working_days[random_index(*rng)];
					for (int p : range(problem.periods.size())) {
						if (problem.periods[p].day == day) {
							dag.disconnect_node(builder.have_worked_node(p));
						}
					}
				}
				consecutive_working_days.clear();
			}
		}
		if (feasible) {
			break;
		}
	}

	if (feasible) {
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
			} else {
				minimum_core_assert(!builder.is_have_worked(node2));
			}
		}
	}

	minimum_core_assert(feasible);
	return feasible;
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
		uniform_real_distribution<double> dual(-10, 100);
		for (int p = 0; p < problem.staff.size(); ++p) {
			for (int r = 0; r < number_of_rows(); ++r) {
				initial_duals[p][r] = dual(rng[p]);
			}
		}
		auto no_fixes =
		    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });

		for (int p = 0; p < problem.staff.size(); ++p) {
			bool success =
			    create_roster_graph(problem, initial_duals[p], p, no_fixes, &solution[p], &rng[p]);
			check(success, "Initial roster generation failed.");
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

		// for (auto var : dual_variables) {
		//	cerr << var << " ";
		//}
		// cerr << endl;

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
