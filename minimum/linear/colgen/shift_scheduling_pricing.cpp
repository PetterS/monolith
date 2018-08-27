#include <fstream>
#include <vector>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/algorithms/dag.h>
#include <minimum/linear/colgen/shift_scheduling_pricing.h>
#include <minimum/linear/proto.h>

namespace minimum {
namespace linear {
namespace colgen {

// Holds a directed graph and keeps track of the mapping from node to (day, shift) and vice
// verca. Can transform the graph into a smaller, equivalent graph.
class ScheduleGraph {
   public:
	ScheduleGraph(const minimum::linear::proto::SchedulingProblem& problem_,
	              int num_nodes_,
	              vector<vector<int>>&& shift_to_node_,
	              vector<int>&& day_off_node_)
	    : dag(num_nodes_),
	      sink(num_nodes_ - 1),
	      problem(problem_),
	      shift_to_node(move(shift_to_node_)),
	      node_to_shift(num_nodes_, -1),
	      node_to_day(num_nodes_, -1),
	      day_off_node_indices(day_off_node_) {
		for (auto d : range(problem.num_days())) {
			for (auto s : range(problem.shift_size())) {
				auto i = node(d, s);
				;
				if (i >= 0) {
					node_to_day[i] = d;
					node_to_shift[i] = s;
				}
			}
		}
	}

	int source_node() const { return source; }
	int sink_node() const { return sink; }
	int node(int day, int shift) const { return shift_to_node[day][shift]; }
	int day_off_node(int day) const { return day_off_node_indices[day]; }
	int shift(int node) const { return node_to_shift[node]; }
	int day(int node) const { return node_to_day[node]; }

	// Given fixed nodes, removes the weights and decreases the minimum and maximum
	// time by the same amount. This makes the resulting problem smaller.
	void optimize_time_limit(const vector<vector<int>>& fixes, int* min, int* max) {
		for (int d : range(problem.num_days())) {
			for (int s : range(problem.shift_size())) {
				if (fixes[d][s] == 1) {
					*min -= problem.shift(s).duration();
					*max -= problem.shift(s).duration();
					auto i = node(d, s);
					minimum_core_assert(i >= 0);
					dag.set_node_weight(i, 0, 0);
				}
			}
		}
	}

	// Modifies the graph in place to a smaller equivalent graph.
	void optimize_size() {
		auto translator = dag.reduce_graph();
		if (translator.made_changes) {
			source = translator.old_to_new[source];
			sink = translator.old_to_new[sink];

			vector<int> new_node_to_shift(dag.size(), -1);
			vector<int> new_node_to_day(dag.size(), -1);
			for (auto i : range(dag.size())) {
				minimum_core_assert(translator.new_to_old[i].size() == 1);
				auto old_i = translator.new_to_old[i][0];
				new_node_to_shift[i] = node_to_shift[old_i];
				new_node_to_day[i] = node_to_day[old_i];
			}
			node_to_shift = move(new_node_to_shift);
			node_to_day = move(new_node_to_day);

			for (int d : range(problem.num_days())) {
				for (auto& i : shift_to_node[d]) {
					if (i >= 0) {
						i = translator.old_to_new[i];
					}
				}
			}
			for (auto& i : day_off_node_indices) {
				if (i >= 0) {
					i = translator.old_to_new[i];
				}
			}
		}
	}

	// Plots the graph in dot format.
	void plot(ostream& out, const vector<vector<int>>& fixes) const {
		auto name = [&](int node) -> string {
			if (node == source) {
				return "Start";
			} else if (node == sink) {
				return "End";
			}

			auto day_off = find(day_off_node_indices.begin(), day_off_node_indices.end(), node);
			if (day_off != day_off_node_indices.end()) {
				return to_string("Off ", day_off - day_off_node_indices.begin());
			}

			auto s = node_to_shift[node];
			return to_string(problem.shift(s).id(), " ", node_to_day[node]);
		};
		auto color = [&](int node) -> string {
			if (node == source) {
				return "goldenrod";
			} else if (node == sink) {
				return "goldenrod";
			}

			auto day_off = find(day_off_node_indices.begin(), day_off_node_indices.end(), node);
			if (day_off != day_off_node_indices.end()) {
				return "yellow";
			}

			auto d = node_to_day[node];
			auto s = node_to_shift[node];
			auto fixed = fixes[d][s];
			if (fixed == 0) {
				// Should not be in the graph.
				return "crimson";
			} else if (fixed == 1) {
				return "chartreuse2";
			}
			return "white";
		};
		auto rank = [&](int node) -> int {
			if (node == source) {
				return -1;
			} else if (node == sink) {
				return problem.num_days();
			}

			auto day_off = find(day_off_node_indices.begin(), day_off_node_indices.end(), node);
			if (day_off != day_off_node_indices.end()) {
				return day_off - day_off_node_indices.begin();
			}

			return node_to_day[node];
		};
		dag.write_dot(out, name, color, rank);
	}

	minimum::algorithms::SortedDAG<double, 2> dag;

   private:
	int source = 0;
	int sink = -1;
	const minimum::linear::proto::SchedulingProblem& problem;
	vector<vector<int>> shift_to_node;
	vector<int> node_to_shift;
	vector<int> node_to_day;
	vector<int> day_off_node_indices;
};

// Builds a ScheduleGraph by keeping track of the number of nodes and their
// associated days and shifts.
class ScheduleGraphBuilder {
   public:
	ScheduleGraphBuilder(const minimum::linear::proto::SchedulingProblem& problem_)
	    : problem(problem_), day_off_node(problem_.num_days()) {
		shift_to_node =
		    make_grid<int>(problem.num_days(), problem.shift_size(), []() { return -1; });
	}

	void add_node(int day, int shift) { shift_to_node[day][shift] = num_nodes++; }
	void add_day_off_node(int day) { day_off_node[day] = num_nodes++; }

	ScheduleGraph build() {
		num_nodes++;  // Sink node.
		return ScheduleGraph(problem, num_nodes, move(shift_to_node), move(day_off_node));
	}

   private:
	int num_nodes = 1;  // Source node only initially.
	const minimum::linear::proto::SchedulingProblem& problem;
	vector<vector<int>> shift_to_node;
	vector<int> day_off_node;
};

bool create_roster_cspp(const minimum::linear::proto::SchedulingProblem& problem,
                        const vector<double>& dual_variables,
                        int p,
                        const vector<vector<int>>& fixes,
                        vector<vector<int>>* solution_for_staff,
                        const std::string& graph_output_file_name) {
	auto& staff = problem.worker(p);

	ScheduleGraphBuilder graph_builder(problem);
	for (auto d : range(problem.num_days())) {
		graph_builder.add_day_off_node(d);

		bool required_day_off =
		    find(staff.required_day_off().begin(), staff.required_day_off().end(), d)
		    != staff.required_day_off().end();
		if (required_day_off) {
			// This is a day off. Do not add any nodes for shifts.
			continue;
		}

		int shift_fixed_to_one = -1;
		for (auto s : range(problem.shift_size())) {
			if (fixes[d][s] == 1) {
				shift_fixed_to_one = s;
				break;
			}
		}
		if (shift_fixed_to_one >= 0) {
			// We only have to add a node for this shift.
			graph_builder.add_node(d, shift_fixed_to_one);
		} else {
			for (auto s : range(problem.shift_size())) {
				if (staff.shift_limit(s).max() >= 1 && fixes[d][s] != 0) {
					graph_builder.add_node(d, s);
				}
			}
		}
	}
	auto graph = graph_builder.build();

	for (auto d : range(problem.num_days())) {
		for (auto s : range(problem.shift_size())) {
			auto i = graph.node(d, s);
			if (i >= 0) {
				// Minutes worked.
				graph.dag.set_node_weight(i, 0, problem.shift(s).duration());
				// This node represents “active” when computing min/max
				// consecutive constraints.
				graph.dag.set_node_weight(i, 1, 1);
			}
		}
	}

	// Cover costs from dual variables.
	int r = problem.worker_size();
	for (const auto& requirement : problem.requirement()) {
		auto i = graph.node(requirement.day(), requirement.shift());
		if (i >= 0) {
			graph.dag.set_node_cost(i, -dual_variables.at(r));
		}
		++r;
	}

	// Shift preferences.
	for (const auto& pref : staff.shift_preference()) {
		auto i = graph.node(pref.day(), pref.shift());
		if (i >= 0) {
			graph.dag.change_node_cost(i, pref.cost());
		}
	}

	// Day-off preferences.
	for (const auto& pref : staff.day_off_preference()) {
		auto i = graph.day_off_node(pref.day());
		if (i >= 0) {
			graph.dag.change_node_cost(i, pref.cost());
		}
	}

	// Connect source and sink.
	graph.dag.add_edge(graph.source_node(), graph.day_off_node(0));
	graph.dag.add_edge(graph.day_off_node(problem.num_days() - 1), graph.sink_node());
	for (auto s : range(problem.shift_size())) {
		int i = graph.node(0, s);
		if (i >= 0) {
			graph.dag.add_edge(graph.source_node(), i);
		}
		i = graph.node(problem.num_days() - 1, s);
		if (i >= 0) {
			graph.dag.add_edge(i, graph.sink_node());
		}
	}

	check(!staff.has_consecutive_days_off_limit()
	          || staff.consecutive_days_off_limit().max() >= problem.num_days(),
	      "A maximum consecutive days off has been set but is not supported.");

	// Connect graph.
	for (auto d : range(problem.num_days() - 1)) {
		// From day off to another day off.
		graph.dag.add_edge(graph.day_off_node(d), graph.day_off_node(d + 1));

		// From shifts to and from day off.
		int min_consecutive = std::max(staff.consecutive_days_off_limit().min(), 1);
		int minimum_day_off = std::min(problem.num_days() - 1, d + min_consecutive);
		bool fixed_working = false;
		for (auto d2 : range(d + 1, minimum_day_off + 1)) {
			for (auto s : range(problem.shift_size())) {
				if (fixes[d2][s] == 1) {
					fixed_working = true;
					break;
				}
			}
		}
		for (auto s : range(problem.shift_size())) {
			int i = graph.node(d + 1, s);
			if (i >= 0) {
				graph.dag.add_edge(graph.day_off_node(d), i);
			}
			i = graph.node(d, s);
			if (!fixed_working && i >= 0) {
				graph.dag.add_edge(i, graph.day_off_node(minimum_day_off));
			}
		}

		// From shifts to other shifts.
		for (auto s1 : range(problem.shift_size())) {
			for (int s2 : problem.shift(s1).ok_after()) {
				int i1 = graph.node(d, s1);
				int i2 = graph.node(d + 1, s2);
				if (i1 >= 0 && i2 >= 0) {
					graph.dag.add_edge(i1, i2);
				}
			}
		}
	}

	// Fixes.
	for (int d : range(problem.num_days())) {
		int one_fixes_for_day = 0;
		for (int s = 0; s < problem.shift_size(); ++s) {
			if (fixes[d][s] == 0) {
				auto i = graph.node(d, s);
				if (i >= 0) {
					graph.dag.disconnect_node(i);
				}
			} else if (fixes[d][s] == 1) {
				one_fixes_for_day++;
				for (int s2 = 0; s2 < problem.shift_size(); ++s2) {
					if (s != s2) {
						auto i = graph.node(d, s2);
						if (i >= 0) {
							graph.dag.disconnect_node(i);
						}
					}
				}
				graph.dag.disconnect_node(graph.day_off_node(d));
			}
		}
		minimum_core_assert(one_fixes_for_day <= 1);
	}

	minimum_core_assert(staff.has_time_limit(),
	                    "This code currently assumes this limit to be set.");
	int time_limit_min = staff.time_limit().min();
	int time_limit_max = staff.time_limit().max();
	graph.optimize_time_limit(fixes, &time_limit_min, &time_limit_max);

	graph.optimize_size();

	if (!graph_output_file_name.empty()) {
		ofstream fout(graph_output_file_name);
		graph.plot(fout, fixes);
	}

	// The shortest path problem does not model everything. Solve it and iterate
	// until feasible.
	vector<int> solution;
	bool is_feasible = false;
	int repetitions = 0;
	do {
		solution.clear();
		minimum_core_assert(staff.has_time_limit() && staff.has_consecutive_shifts_limit(),
		                    "This code currently assumes these limits to be set.");
		minimum::algorithms::resource_constrained_shortest_path(
		    graph.dag,
		    time_limit_min,
		    time_limit_max,
		    staff.consecutive_shifts_limit().min(),
		    staff.consecutive_shifts_limit().max(),
		    &solution);
		minimum_core_assert(!solution.empty());

		for (auto d : range(problem.num_days())) {
			for (auto s : range(problem.shift_size())) {
				(*solution_for_staff)[d][s] = 0;
			}
		}

		for (auto i : solution) {
			auto s = graph.shift(i);
			auto d = graph.day(i);
			if (s >= 0 && d >= 0) {
				(*solution_for_staff)[d][s] = 1;
			}
		}

		is_feasible = true;

		//
		// Maximum shifts.
		//
		vector<pair<int, double>> active_node_costs;
		vector<int> inactive_nodes;
		for (int s = 0; s < problem.shift_size(); ++s) {
			active_node_costs.clear();
			inactive_nodes.clear();
			int working_shifts = 0;
			for (int d : range(problem.num_days())) {
				auto i = graph.node(d, s);
				if (i >= 0) {
					if ((*solution_for_staff)[d][s] == 1) {
						working_shifts += 1;
						if (fixes[d][s] != 1 && graph.dag.get_node(i).cost < 1e9) {
							active_node_costs.emplace_back(i, graph.dag.get_node(i).cost);
						}
					} else {
						// minimum_core_assert(fixes[d][s] != 1);
						inactive_nodes.emplace_back(i);
					}
				}
			}
			if (working_shifts < staff.shift_limit(s).min()) {
				check(
				    false,
				    "Minimum working a specific shift violated but the code does not support that "
				    "constraint.");
			}

			if (working_shifts > staff.shift_limit(s).max()) {
				sort(active_node_costs.begin(), active_node_costs.end(), [](auto& a, auto& b) {
					return a.second > b.second;
				});
				minimum_core_assert(staff.shift_limit(s).max() > 0,
				                    "Disallowed shift should already have been taken care of.");

				is_feasible = false;
				int need_to_drop = working_shifts - staff.shift_limit(s).max();
				// if (repetitions > 5) {
				//	cerr << to_string("Have ", working_shifts, " of shift ", s, ", limit is ",
				// staff.max_shifts.at(s), ". Need to drop ", need_to_drop, " from ",
				// active_node_costs,
				//".\n");
				//}
				// Drop the right amount of shifts to make the solution feasible.
				for (int to_drop = 0; to_drop < need_to_drop; ++to_drop) {
					if (to_drop >= active_node_costs.size()) {
						// Nothing really can be done at this point. Infeasible.
						return false;
					}

					auto i = active_node_costs.at(to_drop).first;
					minimum_core_assert(graph.shift(i) == s, graph.shift(i), " != ", s);

					// Fixes days before and after this day in combination with shifts that cannot
					// follow other shifts can make the problem infeasible if we simply disconnect
					// this node.
					graph.dag.set_node_cost(i, 1e9);
				}
				// Drop all shifts of this type that are not currently used. This can
				// make the process converge faster.
				for (auto i : inactive_nodes) {
					graph.dag.set_node_cost(i, 1e9);
				}
			}
		}

		//
		// Max working weekends.
		//
		struct Day {
			Day() {}
			Day(bool valid_, double node_cost_) : valid(valid_), node_cost(node_cost_) {}
			bool valid = false;
			double node_cost = 0;
		};
		auto working_on_day = [&](int day) -> Day {
			for (int s = 0; s < problem.shift_size(); ++s) {
				if ((*solution_for_staff)[day][s] == 1) {
					auto i = graph.node(day, s);
					return {true, graph.dag.get_node(i).cost};
				}
			}
			return {false, 0};
		};
		auto day_is_fixed_working = [&](int day) -> bool {
			for (int s = 0; s < problem.shift_size(); ++s) {
				if (fixes[day][s] > 0) {
					return true;
				}
			}
			return false;
		};

		struct Weekend {
			Weekend(int saturday_index_, const Day& saturday, const Day& sunday) {
				saturday_index = saturday_index_;
				cost = 0;
				if (saturday.valid) {
					cost += saturday.node_cost;
				}
				if (sunday.valid) {
					cost += sunday.node_cost;
				}
			}
			int saturday_index;
			double cost;
		};

		int working_weekends = 0;
		vector<Weekend> active_weekends;
		vector<Weekend> inactive_weekends;
		for (int d = 5; d < problem.num_days() - 1; d += 7) {
			auto saturday = working_on_day(d);
			auto sunday = working_on_day(d + 1);
			if (saturday.valid || sunday.valid) {
				working_weekends++;
				if (!day_is_fixed_working(d) && !day_is_fixed_working(d + 1)
				    && saturday.node_cost < 1e10 && sunday.node_cost < 1e10) {
					active_weekends.emplace_back(d, saturday, sunday);
				}
			} else {
				inactive_weekends.emplace_back(d, saturday, sunday);
			}
		}

		if (staff.has_working_weekends_limit()
		    && working_weekends < staff.working_weekends_limit().min()) {
			check(
			    false,
			    "Minimum working weekends violated but the code does not support that constraint.");
		}

		if (staff.has_working_weekends_limit()
		    && working_weekends > staff.working_weekends_limit().max()) {
			is_feasible = false;
			sort(active_weekends.begin(),
			     active_weekends.end(),
			     [](const Weekend& a, const Weekend& b) { return a.cost > b.cost; });

			int need_to_drop = working_weekends - staff.working_weekends_limit().max();

			// cerr << to_string(p, " has ", working_weekends, " weekends, limit is ",
			// staff.max_weekends,
			// ". Need to drop ", need_to_drop, " from [ ");
			// for (auto& w : active_weekends) {
			//	cerr << "(" << w.saturday_index << ", " << w.cost << ") ";
			//}
			// cerr << ")\n";

			// Drop the right amount of active weekends.
			for (int to_drop = 0; to_drop < need_to_drop; ++to_drop) {
				if (to_drop >= active_weekends.size()) {
					// Nothing really can be done at this point. Infeasible.
					return false;
				}

				auto& weekend = active_weekends.at(to_drop);
				for (int s = 0; s < problem.shift_size(); ++s) {
					// Fixes days before and after this weekend in combination with shifts that
					// cannot follow other shifts can make the problem infeasible if we simply
					// disconnect this node.
					auto i1 = graph.node(weekend.saturday_index, s);
					if (i1 >= 0) {
						graph.dag.set_node_cost(i1, 1e10);
					}
					auto i2 = graph.node(weekend.saturday_index + 1, s);
					if (i2 >= 0) {
						graph.dag.set_node_cost(i2, 1e10);
					}
				}
			}
			// Forbid all weekends that are currently not used. This can make
			// things converge faster.
			for (const auto& weekend : inactive_weekends) {
				for (int s = 0; s < problem.shift_size(); ++s) {
					auto i1 = graph.node(weekend.saturday_index, s);
					if (i1 >= 0) {
						graph.dag.set_node_cost(i1, 1e10);
					}
					auto i2 = graph.node(weekend.saturday_index + 1, s);
					if (i2 >= 0) {
						graph.dag.set_node_cost(i2, 1e10);
					}
				}
			}
		}

		minimum_core_assert(++repetitions <= 100);
	} while (!is_feasible);

	// Sanity check. Will crash if something has gone very wrong in the graph
	// construction.
	for (auto d : range(problem.num_days())) {
		int num_active = 0;
		for (auto s : range(problem.shift_size())) {
			if ((*solution_for_staff)[d][s] != 0) {
				num_active++;
			}
			if (fixes[d][s] >= 0) {
				minimum_core_assert(fixes[d][s] == (*solution_for_staff)[d][s]);
			}
		}
		minimum_core_assert(num_active <= 1);
	}

	return true;
}
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
