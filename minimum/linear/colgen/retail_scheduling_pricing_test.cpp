#include <vector>

#include <Catch.hpp>

#include <minimum/core/grid.h>
#include <minimum/linear/colgen/retail_scheduling_pricing.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/retail_scheduling.h>

using namespace std;
using namespace minimum::core;
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