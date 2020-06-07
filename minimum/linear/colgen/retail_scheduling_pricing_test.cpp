#include <vector>

#include <catch.hpp>

#include <minimum/core/grid.h>
#include <minimum/core/random.h>
#include <minimum/linear/colgen/retail_scheduling_pricing.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/retail_scheduling.h>

using namespace std;
using namespace minimum::core;
using minimum::linear::colgen::create_retail_shift;
using minimum::linear::colgen::create_roster_graph;
using minimum::linear::colgen::retail_local_search;

minimum::linear::RetailProblem get_problem(std::string name) {
	auto base_dir = minimum::linear::data::get_directory();
	std::ifstream file{base_dir + "/retail/" + name + ".txt"};
	return {file};
}

TEST_CASE("simplest_pricing") {
	auto problem = get_problem("0");
	const auto fixes =
	    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });
	vector<double> duals(problem.staff.size() + problem.periods.size() * problem.num_tasks, 0);
	auto solution = make_grid<int>(problem.periods.size(), problem.num_tasks);
	mt19937 rng;

	for (int staff_index = 0; staff_index < problem.staff.size(); staff_index++) {
		CHECK(create_roster_graph(problem, duals, staff_index, fixes, &solution, &rng));
	}
}

TEST_CASE("simple_pricing") {
	auto problem = get_problem("1_7_10_1");

	const auto fixes =
	    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });
	auto solution = make_grid<int>(problem.periods.size(), problem.num_tasks);
	mt19937 rng;

	for (int p = 0; p < problem.periods.size(); p += 20) {
		vector<double> duals(problem.staff.size() + problem.periods.size(), 0);
		duals[problem.staff.size() + problem.cover_constraint_index(p, 0)] = 1000;
		const int staff_index = 0;
		CHECK(create_roster_graph(problem, duals, staff_index, fixes, &solution, &rng));
		CHECK(solution[p][0] == 1);
		problem.check_feasibility_for_staff(staff_index, solution);
	}
}

TEST_CASE("simple_local_search") {
	auto problem = get_problem("1_7_10_1");
	minimum::linear::colgen::RetailLocalSearchParameters parameters;
	parameters.time_limit_seconds = 0.1;
	retail_local_search(problem, parameters);
}

TEST_CASE("pricing_with_multiple_tasks_and_fixes") {
	auto problem = get_problem("10_7_20_3");

	auto fixes = make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });
	auto solution = make_grid<int>(problem.periods.size(), problem.num_tasks);

	auto rng = seeded_engine<mt19937>();
	uniform_real_distribution<double> eps(-1, 1);
	vector<double> duals(problem.staff.size() + problem.periods.size() * problem.num_tasks, 0);
	for (auto& dual : duals) {
		dual = eps(rng);
	}

	for (int p = 4; p < problem.periods.size() - 4; p += 20) {
		fixes[p][1] = 1;
		fixes[p - 1][2] = 1;

		const int staff_index = 0;
		REQUIRE(create_roster_graph(problem, duals, staff_index, fixes, &solution, &rng));
		CHECK(solution[p][1] == 1);
		CHECK(solution[p - 1][2] == 1);
		problem.check_feasibility_for_staff(staff_index, solution);

		fixes[p][1] = -1;
		fixes[p - 1][2] = -1;
		break;
	}
}

namespace {
void run_existing_solution_test(string problem_name, string solution_name) {
	auto rng = seeded_engine<mt19937>();
	const auto problem = get_problem(problem_name);
	const auto base_dir = minimum::linear::data::get_directory();
	ifstream solution_file{base_dir + "/retail/" + solution_name + ".solution"};
	const auto existing_solution = problem.load_solution(solution_file);
	const auto fixes =
	    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });

	for (int staff_index : range(problem.staff.size())) {
		vector<double> duals(problem.staff.size() + problem.periods.size() * problem.num_tasks, 0);
		auto solution = make_grid<int>(problem.periods.size(), problem.num_tasks);
		for (int period_index : range(problem.periods.size())) {
			for (int t : range(problem.num_tasks)) {
				auto current = existing_solution[staff_index].at(period_index).at(t);
				double dual_change;
				if (current == 1) {
					dual_change = 1000;
				} else {
					dual_change = -1000;
				}
				duals.at(problem.staff.size() + problem.cover_constraint_index(period_index, t)) =
				    dual_change;
			}
		}

		REQUIRE(create_roster_graph(problem, duals, staff_index, fixes, &solution, &rng));

		for (int period_index : range(problem.periods.size())) {
			for (int t : range(problem.num_tasks)) {
				if (solution.at(period_index).at(t)
				    != existing_solution[staff_index][period_index][t]) {
					auto& period = problem.periods[period_index];
					auto& person = problem.staff[staff_index];
					auto msg = to_string("Mismatch for staff ",
					                     person.id,
					                     " at day ",
					                     period.day,
					                     ", ",
					                     period.human_readable_time,
					                     ", task ",
					                     t,
					                     ". Wanted ",
					                     existing_solution[staff_index][period_index][t],
					                     ".");
					FAIL(msg);
				}
			}
		}
	}
}
}  // namespace

TEST_CASE("existing_good_solution_is_feasible_for_pricing") {
	// Very good test that ensures a resonably good solution is actually feasible for the pricer.
	// This small test case has a few good features:
	//   - A single task, so that the heuristic for picking tasks is not tested here.
	//     (should still work though)
	//   - Short enough so that the heuristic for disallowing consecutive days is not used.
	//     (but the duals are created so it should work anyway)
	run_existing_solution_test("1_7_10_1", "test_solution_1a");
}

TEST_CASE("existing_good_solution_with_multiple_shifts_is_feasible_for_pricing") {
	// Same as the previous test but with several tasks.
	run_existing_solution_test("10_7_20_3", "test_solution_10a");
}

TEST_CASE("create_retail_shift_2a") {
	vector<vector<double>> shift_costs = {
	    {0, -1},
	    {0, -1},
	    {0, -1},
	    {0, -1},
	    {0, -1},
	    {-1, 0},
	    {-1, 0},
	    {-1, 0},
	    {-1, 0},
	};
	CHECK(create_retail_shift(shift_costs) == vector<int>{1, 1, 1, 1, 1, 0, 0, 0, 0});
}

TEST_CASE("create_retail_shift_2b") {
	vector<vector<double>> shift_costs = {
	    {0, -10},
	    {0, 1},
	    {0, 1},
	    {0, 1},
	    {0, 1},
	    {0, 1},
	    {0, 1},
	    {0, 1},
	    {0, -1},
	};
	CHECK(create_retail_shift(shift_costs) == vector<int>{1, 1, 1, 1, 0, 0, 0, 0, 0});
}

TEST_CASE("create_retail_shift_2c") {
	vector<vector<double>> shift_costs = {
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	    {10, 0},
	};
	for (int i = 0; i < shift_costs.size(); ++i) {
		auto solution = create_retail_shift(shift_costs);
		CHECK(solution[i] != 0);
		shift_costs[i][0] = -1000;
		solution = create_retail_shift(shift_costs);
		CHECK(solution[i] == 0);
		shift_costs[i][0] = 10;
	}
}

TEST_CASE("create_retail_shift_3a") {
	vector<vector<double>> shift_costs = {
	    {3, 0, 100},
	    {-3, 0, 100},
	    {-3, 0, 100},
	    {-3, 0, 100},
	    {0, 0, 3},
	    {0, 0, -20},
	    {0, 0, 5},
	    {0, 0, 8},
	    {0, 0, 1},
	};
	CHECK(create_retail_shift(shift_costs) == vector<int>{0, 0, 0, 0, 0, 2, 2, 2, 2});
}

TEST_CASE("create_retail_shift_3bc") {
	vector<vector<double>> shift_costs = {
	    {10, 0, 10},
	    {10, -1, 10},
	    {10, -2, 10},
	    {10, 2, 10},
	    {10, -1, 10},
	    {10, -3, 10},
	    {10, 1, 10},
	    {10, 2, 10},
	    {10, -1, 10},
	    {10, -3, 10},
	    {10, 1, 10},
	    {10, -1, 10},
	    {10, -1, 10},
	};
	for (int i = 4; i < shift_costs.size() - 5; ++i) {
		auto solution = create_retail_shift(shift_costs);
		CHECK(solution[i] != 0);
		CHECK(solution[i + 1] != 2);
		shift_costs[i][0] = -1000;
		shift_costs[i + 1][2] = -1000;
		solution = create_retail_shift(shift_costs);
		CHECK(solution[i] == 0);
		CHECK(solution[i + 1] == 2);
		shift_costs[i][0] = 10;
		shift_costs[i + 1][2] = 10;
	}
}
