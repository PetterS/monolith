// Petter Strandmark.
//
// The code in this file requires C++14 generic lambdas.
//

#include <cmath>
#include <limits>
#include <random>

#include <catch.hpp>

#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/dynamic_auto_diff_term.h>
#include <minimum/nonlinear/solver.h>
using namespace minimum::nonlinear;

TEST_CASE("make_term_2") {
	Function function;
	std::vector<double> x = {0, 0};

	auto lambda = [](auto* x) {
		auto d0 = x[1] - x[0] * x[0];
		auto d1 = 1 - x[0];
		return 100 * d0 * d0 + d1 * d1;
	};

	auto term = make_differentiable<2>(lambda);

	function.add_term(term, x.data());

	NewtonSolver solver;
	std::stringstream sout;
	solver.log_function = [&sout](auto str) { sout << str << std::endl; };

	SolverResults results;
	solver.solve(function, &results);
	INFO(sout.str());
	CHECK(std::abs(x[0] - 1.0) < 1e-8);
	CHECK(std::abs(x[1] - 1.0) < 1e-8);
	CHECK(results.exit_success());
}

TEST_CASE("make_term_1_1") {
	Function function;
	double x;
	double y;

	auto lambda = [](auto* x, auto* y) {
		auto d0 = y[0] - x[0] * x[0];
		auto d1 = 1 - x[0];
		return 100 * d0 * d0 + d1 * d1;
	};

	auto term = make_differentiable<1, 1>(lambda);

	function.add_term(term, &x, &y);

	NewtonSolver solver;
	std::stringstream sout;
	solver.log_function = [&sout](auto str) { sout << str << std::endl; };

	SolverResults results;
	solver.solve(function, &results);
	INFO(sout.str());
	CHECK(std::abs(x - 1.0) < 1e-8);
	CHECK(std::abs(y - 1.0) < 1e-8);
	CHECK(results.exit_success());
}
