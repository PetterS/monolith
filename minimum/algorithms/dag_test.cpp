#include <vector>
using namespace std;

#include <catch.hpp>

#include <minimum/algorithms/dag.h>
#include <minimum/core/range.h>
using namespace minimum::algorithms;
using namespace minimum::core;

TEST_CASE("shortest_path_straight") {
	SortedDAG<double, 2> dag(6);
	for (int i : range(dag.size())) {
		dag.set_node_cost(i, i + 1);
		if (i > 0) {
			dag.add_edge(i - 1, i);
		}
	}
	vector<SolutionEntry<double>> solution;
	CHECK(shortest_path(dag, &solution) == 21);

	vector<int> int_solution;
	CHECK(resource_constrained_shortest_path(dag, 0, 10, &int_solution) == 21);
	CHECK(resource_constrained_shortest_path(dag, 0, 10, 2, 6, &int_solution) == 21);

	auto translator = dag.reduce_graph();
	CHECK_FALSE(translator.made_changes);
}

TEST_CASE("shortest_path_small") {
	SortedDAG<double, 2> dag(10);
	for (int i : range(1, 10)) {
		dag.add_edge(i - 1, i, 1.0);
	}
	dag.add_edge(3, 6, 1.0);

	vector<SolutionEntry<double>> solution;
	CHECK(shortest_path(dag, &solution) == 7);
	CHECK(solution[9].prev == 8);
	CHECK(solution[8].prev == 7);
	CHECK(solution[7].prev == 6);
	CHECK(solution[6].prev == 3);
	CHECK(solution[3].prev == 2);
	CHECK(solution[2].prev == 1);
	CHECK(solution[1].prev == 0);
	CHECK(solution[0].prev == -1);

	vector<int> int_solution;
	CHECK(resource_constrained_shortest_path(dag, 0, 5, &int_solution) == 7);
	CHECK(int_solution == (vector<int>{0, 1, 2, 3, 6, 7, 8, 9}));
	CHECK(solution_cost(dag, int_solution) == 7);

	CHECK(resource_constrained_shortest_path(dag, 0, 5, 2, 6, &int_solution) == 7);
	CHECK(int_solution == (vector<int>{0, 1, 2, 3, 6, 7, 8, 9}));
	CHECK(solution_cost(dag, int_solution) == 7);

	auto translator = dag.reduce_graph();
	CAPTURE(to_string(translator.new_to_old));
	CAPTURE(to_string(translator.old_to_new));
	CHECK_FALSE(translator.made_changes);
}

TEST_CASE("shortest_path_0") {
	SortedDAG<double, 2> dag(0);
	vector<SolutionEntry<double>> solution;
	shortest_path(dag, &solution);
	CHECK(solution.empty());

	vector<int> int_solution;
	resource_constrained_shortest_path(dag, 0, 5, &int_solution);
	CHECK(int_solution.empty());
	resource_constrained_shortest_path(dag, 0, 5, 2, 6, &int_solution);
	CHECK(int_solution.empty());
}

TEST_CASE("shortest_path_1") {
	SortedDAG<double, 2> dag(1);
	vector<SolutionEntry<double>> solution;
	CHECK(shortest_path(dag, &solution) == 0);
	CHECK(solution.size() == 1);
	vector<int> int_solution;
	CHECK(resource_constrained_shortest_path(dag, 0, 5, &int_solution) == 0);
	CHECK(int_solution.size() == 1);
	CHECK(resource_constrained_shortest_path(dag, 0, 5, 2, 6, &int_solution) == 0);
	CHECK(int_solution.size() == 1);
}

class ScheduleGraph {
   public:
	int s;
	int t;
	vector<int> day_off_node;
	vector<int> working_node;
	SortedDAG<double, 2> dag;

	ScheduleGraph(int num_days) : dag(1) {
		int num_nodes = 0;
		s = num_nodes++;
		vector<int> tmp_node;
		for (int i : range(num_days)) {
			day_off_node.push_back(num_nodes++);
			working_node.push_back(num_nodes++);
			tmp_node.push_back(num_nodes++);
		}
		t = num_nodes++;
		dag = SortedDAG<double, 2>(num_nodes);

		// Stop nodes for consecutive weights.
		dag.set_node_weight(s, 1, -1);
		dag.set_node_weight(t, 1, -1);

		dag.add_edge(s, day_off_node[0]);
		dag.add_edge(s, working_node[0]);
		dag.add_edge(s, tmp_node[0]);
		for (int i : range(1, num_days)) {
			dag.add_edge(day_off_node[i - 1], tmp_node[i]);
			dag.add_edge(working_node[i - 1], tmp_node[i]);

			dag.add_edge(day_off_node[i - 1], day_off_node[i]);
			dag.add_edge(day_off_node[i - 1], working_node[i]);
			dag.add_edge(working_node[i - 1], day_off_node[i]);
			dag.add_edge(working_node[i - 1], working_node[i]);

			dag.add_edge(tmp_node[i - 1], day_off_node[i]);
			dag.add_edge(tmp_node[i - 1], working_node[i]);
		}

		// Temp node are not meant to be used in the optimal path.
		// But they check that the optimal path is taken.
		for (int i : range(num_days)) {
			dag.set_node_weight(tmp_node[i], 0, 10000);
			dag.set_node_weight(tmp_node[i], 1, 10000);
			dag.set_node_cost(tmp_node[i], 10000);
		}

		dag.add_edge(day_off_node.back(), t);
		dag.add_edge(working_node.back(), t);
		dag.add_edge(tmp_node.back(), t);
	}

	string node_name(int i) const {
		if (i == s) {
			return "s";
		}
		if (i == t) {
			return "t";
		}
		for (auto j : range(day_off_node.size())) {
			if (i == day_off_node[j]) {
				return to_string("off", j);
			}
			if (i == working_node[j]) {
				return to_string("work", j);
			}
		}
		return "err";
	}
};

TEST_CASE("shortest_path_1_resource_constraint") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);

	for (int i : range(num_days)) {
		schedule.dag.set_node_cost(schedule.working_node[i], -1.0);
		schedule.dag.set_node_weight(schedule.working_node[i], 0, 1);
	}
	schedule.dag.set_node_weight(schedule.working_node[4], 0, 2);
	schedule.dag.set_node_weight(schedule.working_node[5], 0, 2);

	vector<vector<SolutionEntry<double>>> solutions;
	CHECK(resource_constrained_shortest_path(schedule.dag, {8}, &solutions) == -8);
	CHECK(solutions.size() <= 3);

	// Dynamic programming.
	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 2, 8, &solution) == -8);
	CHECK(solution_cost(schedule.dag, solution) == -8);
	CHECK(resource_constrained_shortest_path(schedule.dag, 2, 8, 2, 6, &solution) == -8);
	CHECK(solution_cost(schedule.dag, solution) == -8);
}

TEST_CASE("shortest_path_1_resource_constraint_infeasible") {
	SortedDAG<double, 1> dag(10);
	for (int i : range(10)) {
		dag.set_node_weight(i, 0, 1);
		if (i > 0) {
			dag.add_edge(i - 1, i, 1.0);
		}
	}
	vector<int> solution;
	CHECK_THROWS(resource_constrained_shortest_path(dag, 5, 7, &solution));
}

TEST_CASE("shortest_path_max_consecutive") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_cost(schedule.working_node[i], -1.0);
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
	}

	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 0, 0, 2, 6, &solution) == -9);
	CHECK(solution_cost(schedule.dag, solution) == -9);
}

TEST_CASE("shortest_path_min_consecutive") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
	}
	schedule.dag.set_node_cost(schedule.working_node[3], 11.0);
	schedule.dag.set_node_cost(schedule.working_node[4], -100.0);
	schedule.dag.set_node_cost(schedule.working_node[5], -100.0);
	schedule.dag.set_node_cost(schedule.working_node[6], 10.0);

	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 0, 0, 3, 6, &solution) == -190);
	CHECK(solution_cost(schedule.dag, solution) == -190);
}

TEST_CASE("shortest_path_min_max_consecutive") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
	}
	schedule.dag.set_node_cost(schedule.working_node[3], 100);
	schedule.dag.set_node_cost(schedule.working_node[3], 100);
	schedule.dag.set_node_cost(schedule.working_node[4], -100);
	schedule.dag.set_node_cost(schedule.working_node[5], -100);
	schedule.dag.set_node_cost(schedule.working_node[6], 10);
	schedule.dag.set_node_cost(schedule.working_node[7], -1);
	schedule.dag.set_node_cost(schedule.working_node[8], -1);
	schedule.dag.set_node_cost(schedule.working_node[9], 100);

	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 0, 0, 3, 4, &solution) == -191);
	CHECK(solution_cost(schedule.dag, solution) == -191);
}

TEST_CASE("shortest_path_resource_and_consecutive") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_cost(schedule.working_node[i], -1.0);
		schedule.dag.set_node_weight(schedule.working_node[i], 0, 1);
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
	}
	schedule.dag.set_node_cost(schedule.working_node[3], 10.0);
	schedule.dag.set_node_cost(schedule.working_node[4], -100.0);
	schedule.dag.set_node_cost(schedule.working_node[5], -100.0);
	schedule.dag.set_node_cost(schedule.working_node[6], 10.0);

	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 6, 6, &solution) == -204);
	CHECK(solution_cost(schedule.dag, solution) == -204);
	CHECK(resource_constrained_shortest_path(schedule.dag, 6, 6, 3, 10, &solution) == -193);
	CHECK(solution_cost(schedule.dag, solution) == -193);

	auto translator = schedule.dag.reduce_graph();
	CAPTURE(to_string(translator.new_to_old));
	CAPTURE(to_string(translator.old_to_new));
	CHECK_FALSE(translator.made_changes);
}

TEST_CASE("shortest_path_min_consecutive_endpoints") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
		schedule.dag.set_node_cost(schedule.working_node[i], 100);
	}
	schedule.dag.set_node_cost(schedule.working_node[0], -1);
	schedule.dag.set_node_cost(schedule.working_node[num_days - 1], -1);

	vector<int> solution;
	// TODO: Allow less then minimum consecutive in the beginning.
	// CHECK(resource_constrained_shortest_path(schedule.dag, 0, 0, 2, 5, &solution) == -2);
	// CHECK(solution_cost(schedule.dag, solution) == -2);
}

TEST_CASE("shortest_path_max_consecutive_endpoints") {
	int num_days = 10;
	ScheduleGraph schedule(num_days);
	for (int i : range(num_days)) {
		schedule.dag.set_node_weight(schedule.working_node[i], 1, 1);
		schedule.dag.set_node_cost(schedule.working_node[i], 100);
	}
	schedule.dag.set_node_cost(schedule.working_node[0], -1);
	schedule.dag.set_node_cost(schedule.working_node[1], -1);
	schedule.dag.set_node_cost(schedule.working_node[2], -1);
	schedule.dag.set_node_cost(schedule.working_node[num_days - 3], -1);
	schedule.dag.set_node_cost(schedule.working_node[num_days - 2], -1);
	schedule.dag.set_node_cost(schedule.working_node[num_days - 1], -1);

	vector<int> solution;
	CHECK(resource_constrained_shortest_path(schedule.dag, 0, 0, 2, 3, &solution) == -6);
	CHECK(solution_cost(schedule.dag, solution) == -6);
}

TEST_CASE("subproblems") {
	int num_days = 20;
	ScheduleGraph schedule(num_days);

	for (int i : range(num_days)) {
		schedule.dag.set_node_cost(schedule.working_node[i], -1.0);
		schedule.dag.set_node_weight(schedule.working_node[i], 0, 1);
	}
	schedule.dag.set_node_weight(schedule.working_node[4], 0, 2);
	schedule.dag.set_node_weight(schedule.working_node[5], 0, 2);
	schedule.dag.set_node_weight(schedule.working_node[14], 0, 2);
	schedule.dag.set_node_weight(schedule.working_node[15], 0, 2);

	vector<int> solution;
	resource_constrained_shortest_path_partial(schedule.dag, 2, 8, 2, 6, 2, &solution);
	CHECK(solution_cost(schedule.dag, solution) == -8);

	resource_constrained_shortest_path_partial(schedule.dag, 2, 8, 2, 2, 2, &solution);
}

TEST_CASE("reduce_graph") {
	SortedDAG<double, 1> dag(10);
	dag.add_edge(0, 5);
	dag.add_edge(5, 9);

	dag.add_edge(5, 6);
	dag.add_edge(6, 7);
	dag.add_edge(6, 8);
	dag.add_edge(1, 2);
	dag.add_edge(2, 3);

	auto translator = dag.reduce_graph();
	CAPTURE(to_string(translator.new_to_old));
	CAPTURE(to_string(translator.old_to_new));

	REQUIRE(dag.size() == 3);

	REQUIRE(translator.old_to_new.size() == 10);
	CHECK(translator.old_to_new.at(0) == 0);
	CHECK(translator.old_to_new.at(1) == -1);
	CHECK(translator.old_to_new.at(5) == 1);
	CHECK(translator.old_to_new.at(9) == 2);

	REQUIRE(translator.new_to_old.size() == 3);
	CHECK(translator.new_to_old.at(0) == vector<int>{0});
	CHECK(translator.new_to_old.at(1) == vector<int>{5});
	CHECK(translator.new_to_old.at(2) == vector<int>{9});
}

TEST_CASE("merge_graph") {
	SortedDAG<double, 1> dag(17);

	// {0, 1, 2, 3} expected to be merged into one node.
	dag.add_edge(0, 1);
	dag.add_edge(1, 2);
	dag.add_edge(2, 3);

	dag.add_edge(3, 4);
	dag.add_edge(3, 5);
	dag.add_edge(4, 5);
	dag.add_edge(4, 6);
	dag.add_edge(5, 6);

	// {6, 7, 8, 9} expected to be merged into one node.
	dag.add_edge(6, 7);
	dag.add_edge(7, 8);
	dag.add_edge(8, 9);

	dag.add_edge(9, 10);
	dag.add_edge(9, 11);
	dag.add_edge(10, 11);
	dag.add_edge(10, 12);
	dag.add_edge(11, 12);

	// {12, 13, 14, 15, 16} expected to be merged into one node.
	dag.add_edge(12, 13);
	dag.add_edge(13, 14);
	dag.add_edge(14, 15);
	dag.add_edge(15, 16);

	auto translator = dag.merge_graph([](int) { return true; });

	REQUIRE(dag.size() == 17 - 3 - 3 - 4);

	REQUIRE(translator.old_to_new.size() == 17);
	CHECK(translator.old_to_new[0] == translator.old_to_new[1]);
	CHECK(translator.old_to_new[0] == translator.old_to_new[2]);
	CHECK(translator.old_to_new[0] == translator.old_to_new[3]);

	CHECK(translator.old_to_new[6] == translator.old_to_new[7]);
	CHECK(translator.old_to_new[6] == translator.old_to_new[8]);
	CHECK(translator.old_to_new[6] == translator.old_to_new[9]);

	CHECK(translator.old_to_new[12] == translator.old_to_new[13]);
	CHECK(translator.old_to_new[12] == translator.old_to_new[14]);
	CHECK(translator.old_to_new[12] == translator.old_to_new[15]);
	CHECK(translator.old_to_new[12] == translator.old_to_new[16]);

	set<int> indices = {translator.old_to_new[0],
	                    translator.old_to_new[4],
	                    translator.old_to_new[5],
	                    translator.old_to_new[6],
	                    translator.old_to_new[10],
	                    translator.old_to_new[11],
	                    translator.old_to_new[12]};
	CHECK(indices.size() == dag.size());

	REQUIRE(translator.new_to_old.size() == dag.size());

	CHECK((translator.new_to_old[0] == vector<int>{0, 1, 2, 3}));
	CHECK((translator.new_to_old[1] == vector<int>{4}));
	CHECK((translator.new_to_old[2] == vector<int>{5}));
	CHECK((translator.new_to_old[3] == vector<int>{6, 7, 8, 9}));
	CHECK((translator.new_to_old[4] == vector<int>{10}));
	CHECK((translator.new_to_old[5] == vector<int>{11}));
	CHECK((translator.new_to_old[6] == vector<int>{12, 13, 14, 15, 16}));

	CHECK(dag.get_node(0).edges.size() == 2);
}

TEST_CASE("shortest_path_1_resource_constraint_with_edge_weights") {
	const int n = 10;
	SortedDAG<int, 1, 1> dag(n);

	//  0 --> 2 --> 4 --> 6 --> 8
	//  V  /  V  /  V  /  V  /  V
	//  1 --> 3 --> 5 --> 7 --> 9

	// All nodes cost -1, and all edges have weight 1.

	for (int i : range(n)) {
		if (i + 1 < n) {
			auto& edge = dag.add_edge(i, i + 1, 0);
			edge.weights[0] = 1;
		}
		if (i + 2 < n) {
			auto& edge = dag.add_edge(i, i + 2, 0);
			edge.weights[0] = 1;
		}
		dag.set_node_cost(i, -1);
	}

	std::vector<int> solution;
	auto value = resource_constrained_shortest_path(dag, 0, 10, &solution);
	CHECK(value == -n);
	value = resource_constrained_shortest_path(dag, 0, 7, &solution);
	CHECK(value == -8);
	CHECK_THROWS(resource_constrained_shortest_path(dag, 0, 3, &solution));

	// Now all nodes instead cost 1.
	for (int i : range(n)) {
		dag.set_node_cost(i, 1);
	}

	value = resource_constrained_shortest_path(dag, 0, 10, &solution);
	CHECK(value == 6);
	value = resource_constrained_shortest_path(dag, 8, 10, &solution);
	CHECK(value == 9);
	CHECK_THROWS(resource_constrained_shortest_path(dag, 11, 14, &solution));
}
