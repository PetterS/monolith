#include <list>
#include <map>

#include <catch.hpp>

#include <minimum/core/numeric.h>

using namespace minimum::core;

TEST_CASE("is_feasible") {
	CHECK(is_feasible(0., 0.5, 1., 1e-6));
	CHECK(is_feasible(0., 0.5, 1., 0.));

	CHECK(is_feasible(0., 1000 + 1e-3, 1000., 1e-6));
	CHECK(is_feasible(0., 1 + 1e-6, 1., 1e-6));
	CHECK(is_feasible(0., -1e-6, 1., 1e-6));

	CHECK_FALSE(is_feasible(0., 1000 + 1e-3, 1000., 1e-7));
	CHECK_FALSE(is_feasible(0., 1 + 1e-6, 1., 1e-7));
	CHECK_FALSE(is_feasible(0., -1e-6, 1., 1e-7));
}

TEST_CASE("is_equal") {
	CHECK(is_equal(0., 0., 1e-6));
	CHECK(is_equal(0., 0., 0.));

	CHECK(is_equal(1e-6, 0., 1e-6));
	CHECK(is_equal(1e-6, 0., 1e-6));
	CHECK(is_equal(1 + 1e-6, 1., 1e-6));
	CHECK(is_equal(1 - 1e-6, 1., 1e-6));
	CHECK(is_equal(1000 + 1e-3, 1000., 1e-6));
	CHECK(is_equal(1000 - 1e-3, 1000., 1e-6));

	CHECK_FALSE(is_equal(1e-6, 0., 1e-7));
	CHECK_FALSE(is_equal(1e-6, 0., 1e-7));
	CHECK_FALSE(is_equal(1 + 1e-6, 1., 1e-7));
	CHECK_FALSE(is_equal(1 - 1e-6, 1., 1e-7));
	CHECK_FALSE(is_equal(1000 + 1e-3, 1000., 1e-7));
	CHECK_FALSE(is_equal(1000 - 1e-3, 1000., 1e-7));
}

TEST_CASE("integer_residual") {
	for (auto sign : {1.0, -1.0}) {
		CHECK(integer_residual(sign * 0.5) == Approx(0.5));
		CHECK(integer_residual(sign * 0.9) == Approx(0.1));
		CHECK(integer_residual(sign * 10.0) == 0);
	}
}
