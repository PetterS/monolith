#include <random>
#include <vector>

#include <minimum/algorithms/dag.h>
#include <minimum/core/numeric.h>
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

		vector<double> best_dual(problem.periods.size());
		for (auto p : range(problem.periods.size())) {
			best_dual[p] = numeric_limits<double>::min();
			for (auto t : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(p, t);
				int row = problem.staff.size() + c;
				double dual = duals.at(row);
				if (dual > best_dual[p]) {
					best_dual[p] = dual;
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
		int prev_fix1 = -1;
		int prev_fix2 = -1;
		int prev_fix3 = -1;
		int prev_fix4 = -1;
		for (int p : range(problem.periods.size())) {
			int this_t = -1;
			for (auto t : range(problem.num_tasks)) {
				if (fixes[p][t] >= 0) {
					minimum_core_assert(fixes[p][t] == 1, "Can not handle zero fixes yet.");
					minimum_core_assert(any_fix[p] == 0);
					any_fix[p] = 1;
					this_t = t;
				}
			}

			// This is subtle!
			//
			// If  period p   is fixed to T1
			// and period p-1 is fixed to T2
			//
			// We need to count p-1, p-2, p-3 and p-4 as fixes as well.
			// Otherwise there might not be room to cram in a full one-hour shift of T2.
			if (this_t >= 0
			    && ((prev_fix4 >= 0 && this_t != prev_fix4)
			        || (prev_fix3 >= 0 && this_t != prev_fix3)
			        || (prev_fix2 >= 0 && this_t != prev_fix2)
			        || (prev_fix1 >= 0 && this_t != prev_fix1))) {
				any_fix.at(p - 4) = 1;
				any_fix.at(p - 3) = 1;
				any_fix.at(p - 2) = 1;
				any_fix.at(p - 1) = 1;
			}
			prev_fix4 = prev_fix3;
			prev_fix3 = prev_fix2;
			prev_fix2 = prev_fix1;
			prev_fix1 = this_t;
			// End subtle part.
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

// Given an array of costs for all tasks in a shift, create the optimal shift.
// The minimum length of 4 is hard-coded in this function.
vector<int> create_retail_shift(const vector<vector<double>>& shift_costs) {
	const int length = shift_costs.size();
	const int num_tasks = shift_costs.at(0).size();
	const int source = 0;
	const int sink = length * num_tasks + 1;
	minimum_core_assert(length >= 5);
	for (int i : range(length)) {
		for (int t : range(num_tasks)) {
			auto val = shift_costs.at(i).at(t);
			minimum_core_assert(-1e10 <= val && val <= 1e10);
		}
	}

	minimum::algorithms::SortedDAG<double, 0, 0> dag(length * num_tasks + 2);

	auto node = [num_tasks](int i, int task) { return i * num_tasks + task + 1; };

	for (int t : range(num_tasks)) {
		dag.add_edge(source,
		             node(3, t),
		             shift_costs.at(0).at(t) + shift_costs.at(1).at(t) + shift_costs.at(2).at(t)
		                 + shift_costs.at(3).at(t));
		dag.add_edge(node(length - 1, t), sink);
	}

	for (int i : range(length - 1)) {
		for (int t : range(num_tasks)) {
			// Can always stay on the same task.
			dag.add_edge(node(i, t), node(i + 1, t), shift_costs.at(i + 1).at(t));
			// Switching tasks require waiting 1h.
			if (i < length - 4) {
				for (int t2 : range(num_tasks)) {
					if (t != t2) {
						double cost = shift_costs.at(i + 1).at(t2) + shift_costs.at(i + 2).at(t2)
						              + shift_costs.at(i + 3).at(t2) + shift_costs.at(i + 4).at(t2);
						dag.add_edge(node(i, t), node(i + 4, t2), cost);
					}
				}
			}
		}
	}

	vector<minimum::algorithms::SolutionEntry<double>> raw_solution;
	auto value = minimum::algorithms::shortest_path(dag, &raw_solution);
	vector<int> solution(length, -1);
	int node_index = raw_solution.size() - 1;
	// Skip sink.
	node_index = raw_solution[node_index].prev;
	int i_prev = -1;
	int t_prev = -1;
	while (node_index >= 0) {
		int i = -1;
		int t = -1;
		if (source < node_index) {
			i = (node_index - 1) / num_tasks;
			t = (node_index - 1) % num_tasks;
		}
		if (t_prev >= 0) {
			for (int idx = i + 1; idx <= i_prev; ++idx) {
				solution[idx] = t_prev;
			}
		}
		i_prev = i;
		t_prev = t;

		node_index = raw_solution[node_index].prev;
	}

	double computed_value = 0;
	for (int i : range(length)) {
		auto t = solution[i];
		minimum_core_assert(0 <= t && t < num_tasks);
		computed_value += shift_costs[i][t];
	}
	minimum_core_assert(relative_error(computed_value, value) < 1e-7);
	// cerr << "-- Value is " << value << endl;
	return solution;
}

void remove_days_with_fixes(const RetailProblem& problem,
                            vector<int>* days,
                            vector<vector<int>> fixes) {
	for (auto itr = days->begin(); itr != days->end(); ++itr) {
		bool day_ok = true;
		for (int p = problem.day_start(*itr); p < problem.day_start(*itr + 1) && day_ok; ++p) {
			for (auto t : range(problem.num_tasks)) {
				if (fixes[p][t] == 1) {
					day_ok = false;
					break;
				}
			}
		}
		if (!day_ok) {
			itr = days->erase(itr);
			itr--;
		}
	}
}

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
		vector<int> consecutive_working_days;
		consecutive_working_days.reserve(problem.num_days);
		for (int day : range(problem.num_days + 1)) {
			if (working_on_day[day] == 1) {
				consecutive_working_days.push_back(day);
			} else {
				if (consecutive_working_days.size() > 5) {
					feasible = false;

					remove_days_with_fixes(problem, &consecutive_working_days, fixes);
					uniform_int_distribution<int> random_index(0,
					                                           consecutive_working_days.size() - 1);
					int day_to_remove = consecutive_working_days[random_index(*rng)];
					for (int p : range(problem.periods.size())) {
						if (problem.periods[p].day == day_to_remove) {
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
		if (problem.num_tasks == 1) {
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
		} else {
			for (int i = 0; i + 1 < graph_solution.size(); ++i) {
				const int node1 = graph_solution[i];
				const int node2 = graph_solution[i + 1];
				if (builder.is_not_worked(node1) && builder.is_have_worked(node2)) {
					const int p1 = builder.get_period(node1);
					const int p2 = builder.get_period(node2);

					vector<vector<double>> shift_costs(p2 - p1 + 1,
					                                   vector<double>(problem.num_tasks, 0.0));
					for (int p = p1; p <= p2; ++p) {
						bool has_fix = false;
						for (int t : range(problem.num_tasks)) {
							if (fixes.at(p).at(t) == 1) {
								shift_costs[p - p1][t] = -1e7;
								minimum_core_assert(!has_fix);
								has_fix = true;
							} else if (fixes.at(p).at(t) == 0) {
								shift_costs[p - p1][t] = 1e7;
								minimum_core_assert(!has_fix);
								has_fix = true;
							} else {
								int c = problem.cover_constraint_index(p, t);
								int row = problem.staff.size() + c;
								double dual = duals.at(row);
								shift_costs[p - p1][t] = -dual;
							}
						}
					}

					const auto shift_solution = create_retail_shift(shift_costs);
					minimum_core_assert(shift_solution.size() == p2 - p1 + 1);
					for (int p = p1, sol_index = 0; p <= p2; ++p, ++sol_index) {
						const int picked_task = shift_solution[sol_index];
						(*solution)[p][picked_task] = 1;
						for (int t : range(problem.num_tasks)) {
							if (fixes.at(p).at(t) == 1) {
								if (picked_task != t) {
									// The solver did not respect our fixes? This can happen in the
									// following situation. Consider the following fixes:
									//
									// time →
									//
									// *  1  1  1  *  0  0  0  1  ...
									//    ↑
									//    |
									//    | Start of shift
									//
									// (from a real case)
									//
									// There is no way to respect these fixes with the given shift
									// start. A more common scenario is prevented by the "subtle"
									// code at the start of this function.
									//
									// We could use a SAT solver to prove that the second star has
									// to be a 0 before solving the first graph problem, but let's
									// leave that for future work.
									//
									// We can not create a legal shift here so we abort.
									cerr << "-- Pricing: unable to handle fixes in period " << p
									     << " task " << t << " picked " << picked_task
									     << ". Staff is " << staff_index << ".\n";
									return false;
								}
							} else if (fixes.at(p).at(t) == 0) {
								minimum_core_assert(picked_task != t);
							}
						}
					}
				} else {
					minimum_core_assert(!builder.is_have_worked(node2));
				}
			}
		}

		problem.check_feasibility_for_staff(staff_index, *solution);
	}

	return feasible;
}

}  // namespace colgen
}  // namespace linear
}  // namespace minimum