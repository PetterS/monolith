#include <fstream>

#include <catch.hpp>

#include <minimum/linear/data/util.h>
#include <minimum/linear/retail_scheduling.h>

minimum::linear::RetailProblem get_problem(std::string name) {
	auto base_dir = minimum::linear::data::get_directory();
	std::ifstream file{base_dir + "/retail/" + name + ".txt"};
	return {file};
}

TEST_CASE("1_7_10_1") {
	auto problem = get_problem("1_7_10_1");
	CHECK(problem.num_days == 7);
}

TEST_CASE("107_14_70_4") {
	auto problem = get_problem("107_14_70_4");
	CHECK(problem.num_days == 14);
	CHECK(problem.staff.size() == 70);
	CHECK(problem.num_tasks == 4);
}
