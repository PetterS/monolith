#include <fstream>

#include <catch.hpp>

#include <minimum/core/range.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/retail_scheduling.h>

using namespace minimum::core;

minimum::linear::RetailProblem get_problem(std::string name) {
	auto base_dir = minimum::linear::data::get_directory();
	return {std::ifstream{base_dir + "/retail/" + name + ".txt"}};
}

TEST_CASE("1_7_10_1") {
	auto problem = get_problem("1_7_10_1");
	CHECK(problem.num_days == 7);
}

TEST_CASE("107_14_70_4") {
	const auto problem = get_problem("107_14_70_4");
	CHECK(problem.num_days == 14);
	CHECK(problem.staff.size() == 70);
	CHECK(problem.num_tasks == 4);

	for (int c : range(problem.num_cover_constraints())) {
		int d = problem.cover_constraint_to_period(c);
		int s = problem.cover_constraint_to_task(c);
		REQUIRE(problem.cover_constraint_index(d, s) == c);
	}
}
