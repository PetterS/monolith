#include <random>
#include <vector>

#include <minimum/algorithms/dag.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/retail_scheduling_pricing.h>

using namespace std;
using namespace minimum::core;

namespace minimum {
namespace linear {
namespace colgen {

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

		// The subset sum structure is used to calculate sums of dual variables for shifts.
		vector<double> subset_sum(problem.periods.size());
		subset_sum[0] = best_dual[0];
		for (auto p : range(size_t{1}, problem.periods.size())) {
			subset_sum[p] = subset_sum[p - 1] + best_dual[p];
		}

		// Data structures respresenting fixes.
		vector<int> any_fix(problem.periods.size(), 0);
		for (int p : range(problem.periods.size())) {
			for (auto t : range(problem.num_tasks)) {
				if (fixes[p][t] >= 0) {
					minimum_core_assert(fixes[p][t] == 1, "Can not handle zero fixes yet.");
					any_fix[p] = 1;
				}
			}
		}
		vector<int> any_fix_subset_sum(problem.periods.size());
		any_fix_subset_sum[0] = any_fix[0];
		for (auto p : range(size_t{1}, problem.periods.size())) {
			any_fix_subset_sum[p] = any_fix_subset_sum[p - 1] + any_fix[p];
		}

		minimum::algorithms::SortedDAG<double, 1, 1> dag(num_nodes());

		dag.add_edge(source(), not_worked_node(0));
		for (int p : range(problem.periods.size() - 1)) {
			if (any_fix[p] != 1) {
				dag.add_edge(not_worked_node(p), not_worked_node(p + 1));
			}
		}
		if (any_fix[problem.periods.size() - 1] != 1) {
			dag.add_edge(not_worked_node(problem.periods.size() - 1), sink());
		}

		for (int day = 0; day < problem.num_days; ++day) {
			int day_end = problem.day_start(day + 1);
			for (int p = problem.day_start(day); p < day_end; ++p) {
				auto& period = problem.periods.at(p);

				if (period.can_start) {
					for (int shift_length = 6 * 4; shift_length <= 10 * 4; ++shift_length) {
						int last_node_in_shift = p + shift_length - 1;
						if (last_node_in_shift >= problem.periods.size()) {
							break;
						}

						double shift_value = subset_sum[last_node_in_shift];
						if (p > 0) {
							shift_value -= subset_sum[p - 1];
						}

						auto& edge = dag.add_edge(
						    not_worked_node(p), have_worked_node(last_node_in_shift), -shift_value);
						edge.weights[0] = shift_length;
					}
				}

				int next_day_in_14_hours = p + 14 * 4 + 1;
				if (next_day_in_14_hours < day_end) {
					next_day_in_14_hours = day_end;
				}
				int fixes_in_period = any_fix_subset_sum.at(min<int>(next_day_in_14_hours - 1,
				                                                     problem.periods.size() - 1))
				                      - any_fix_subset_sum[p];
				minimum_core_assert(fixes_in_period >= 0);

				if (fixes_in_period == 0) {
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
	for (int attempts = 1; attempts <= 10; ++attempts) {
		try {
			auto cost =
			    resource_constrained_shortest_path(dag,
			                                       problem.staff.at(staff_index).min_minutes / 15,
			                                       problem.staff.at(staff_index).max_minutes / 15,
			                                       &graph_solution);
		} catch (std::exception& e) {
			if (attempts == 1) {
				throw std::runtime_error(
				    to_string("Shortest path failed in attempt ", attempts, ": ", e.what()));
			} else {
				// This can happen when we forbid a day. TODO: Forbid another day when the picked
				// day has a fix.
				return false;
			}
		}
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

}  // namespace colgen
}  // namespace linear
}  // namespace minimum