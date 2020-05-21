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

minimum::linear::RetailProblem get_problem(std::string name) {
	auto base_dir = minimum::linear::data::get_directory();
	std::ifstream file{base_dir + "/retail/" + name + ".txt"};
	return {file};
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
