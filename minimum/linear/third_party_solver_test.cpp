#include <catch.hpp>

#include <minimum/linear/glpk.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/testing.h>
using namespace minimum::linear;

TEST_CASE("glpk_ip") {
	GlpkSolver solver;
	solver.set_silent(true);
	CHECK_NOTHROW(simple_integer_programming_tests(solver));
}

TEST_CASE("glpk_lp") {
	GlpkSolver solver;
	solver.set_silent(true);
	CHECK_NOTHROW(simple_linear_programming_tests(solver));
}

TEST_CASE("scs") {
	ScsSolver solver;
	solver.set_convergence_tolerance(1e-9);
	solver.set_silent(true);
	CHECK_NOTHROW(simple_linear_programming_tests(solver));
}

TEST_CASE("scs-warm-start") {
	// TODO: Fix SCS warm start.
	if (false) {
		IP ip;
		create_soduku_IP(ip, 3);
		ScsSolver solver;
		solver.set_convergence_tolerance(1e-9);
		// solver.set_silent(true);
		solver.set_max_iterations(10'000);

		REQUIRE(solver.solutions(&ip)->get());
		REQUIRE(ip.is_feasible());

		solver.set_max_iterations(2);
		CHECK(solver.solutions(&ip)->get());
		CHECK(ip.is_feasible());
	}
}
