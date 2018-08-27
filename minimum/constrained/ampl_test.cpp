#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <streambuf>
#include <string>

#include <catch.hpp>

#include <minimum/constrained/ampl.h>
#include <minimum/constrained/successive.h>
#include <minimum/core/check.h>
#include <minimum/core/numeric.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/data/util.h>
#include <minimum/nonlinear/solver.h>
using minimum::core::check;
using minimum::core::range;
using minimum::nonlinear::data::get_directory;

namespace {
std::string get_file(const std::string& filename) {
	std::ifstream file(filename, std::ios::binary);
	check(!!file, "Could not open ", filename);
	return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void test_problem(std::string_view data, double expected_objective) {
	minimum::constrained::AmplParser parser(data);
	auto problem = parser.parse_ampl();

	if (!problem.has_constraints()) {
		minimum::nonlinear::Function f;
		problem.create_function(&f);
		return;
	}

	minimum::nonlinear::ConstrainedFunction f;
	problem.create_function(&f);

	std::stringstream constraint_log;
	for (auto& entry : f.constraints()) {
		constraint_log << entry.first << "\n";
	}
	{
		INFO(constraint_log.str());
		CHECK(!f.constraints().empty());
	}

	minimum::constrained::SuccessiveLinearProgrammingSolver solver;
	solver.function_improvement_tolerance = 1e-6;
	std::stringstream log;
	solver.log_function = [&log](std::string_view s) { log << s << std::endl; };

	auto result = solver.solve(&f);
	for (auto i : range(problem.get_variables().size())) {
		log << problem.get_variables()[i].name << "=" << problem.get_values()[i] << "\n";
	}
	{
		INFO(log.str());
		CHECK(result == minimum::constrained::SuccessiveLinearProgrammingSolver::Result::OK);
	}
	CHECK(f.is_feasible());
	if (!std::isnan(expected_objective)) {
		auto objective = f.objective().evaluate();
		CAPTURE(expected_objective);
		CAPTURE(objective);
		CHECK(minimum::core::relative_error(expected_objective, objective) <= 1e-4);
	}
}

void test_global_problem(std::string problem_name,
                         double expected_objective = std::numeric_limits<double>::quiet_NaN()) {
	std::string data = get_file(get_directory() + "/coconut/global/" + problem_name + ".mod");
	test_problem(data, expected_objective);
}
}  // namespace

TEST_CASE("circle") {
	test_problem(R"(
		var x := 3;
		var y := 4;
		minimize obj: 0;
		subject to
			first: x^2 + y^2 <= 1;
			second: (x - 2)^2 + y^2 <= 1;
	)",
	             0.0);
}

TEST_CASE("ex14_1_3") { test_global_problem("ex14_1_3", 0.0); }
TEST_CASE("abel") { test_global_problem("abel", 225.19460000); }
TEST_CASE("alkyl") { test_global_problem("alkyl", -1.76500000); }
TEST_CASE("bearing") { test_global_problem("bearing"); }

// Unconstrained problem. Only test parsing + creation.
TEST_CASE("ex4_1_5") { test_global_problem("ex4_1_5"); }