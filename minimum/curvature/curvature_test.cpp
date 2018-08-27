// Petter Strandmark 2013.

#include <cmath>
#include <random>

#include <badiff.h>
#include <fadiff.h>

#include <catch.hpp>

#include <minimum/curvature/curvature.h>

using namespace minimum::curvature;

#define EXPECT_NEAR(a, b, tol)                                                      \
	if (std::abs(a) >= tol || std::abs(b) >= tol) {                                 \
		CHECK((std::abs((a) - (b)) / std::max(std::abs(a), std::abs(b))) <= (tol)); \
	}

TEST_CASE("compute_curvature/helix", "") {
	double x1 = 1;
	double y1 = 0;
	double z1 = 0;

	double x2 = 0.885456025653210;
	double y2 = 0.464723172043769;
	double z2 = 0.025641025641026;

	double x3 = 0.568064746731156;
	double y3 = 0.822983865893656;
	double z3 = 0.051282051282051;

	// Real value   0.481288653714417
	// Approx value 0.495821512020759;

	const double power = 2.0;
	const int n_points = 1000000;
	double k2_int_pair =
	    compute_curvature(x1, y1, z1, x2, y2, z2, x3, y3, z3, power, false, n_points);
	CHECK(curvature_cache_misses == 1);
	EXPECT_NEAR(k2_int_pair, 0.495821512020759, 1e-6);

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		k2_int_pair = compute_curvature(x1, y1, z1, x2, y2, z2, x3, y3, z3, power, true, n_points);
		CHECK(curvature_cache_misses == 2);
		CHECK(curvature_cache_hits == iter - 1);
	}
}

TEST_CASE("compute_torsion/helix", "") {
	double x1 = 1;
	double y1 = 0;
	double z1 = 0;

	double x2 = 0.885456025653210;
	double y2 = 0.464723172043769;
	double z2 = 0.025641025641026;

	double x3 = 0.568064746731156;
	double y3 = 0.822983865893656;
	double z3 = 0.051282051282051;

	double x4 = 0.120536680255323;
	double y4 = 0.992708874098054;
	double z4 = 0.076923076923077;

	const int n_points = 1000000;

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		double t1_int_pair =
		    compute_torsion(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 1.0, true, n_points);
		EXPECT_NEAR(t1_int_pair, 0.026620592104080, 1e-6);
	}

	// Iterate to test cache.
	for (int iter = 1; iter <= 10; ++iter) {
		double t2_int_pair =
		    compute_torsion(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, 2.0, true, n_points);
		EXPECT_NEAR(t2_int_pair, 0.001522684996624, 1e-6);
	}
}

TEST_CASE("compute_torsion/crooked_line_in_plane", "") {
	float test = compute_torsion<float>(3.0, 4.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0, 1.0, 4.0, 2.0, 0.0);
	CHECK(test == test);
	EXPECT_NEAR(test, 0.0f, 1e-6f);
}

TEST_CASE("differentiate_curvature") {
	typedef fadbad::F<double, 9> Dual;
	double h = 1e-3;

	Dual x1 = 1;
	x1.diff(0);
	Dual y1 = 0;
	y1.diff(1);
	Dual z1 = 0;
	z1.diff(2);

	Dual x2 = 0.885456025653210;
	x2.diff(3);
	Dual y2 = 0.464723172043769;
	y2.diff(4);
	Dual z2 = 0.025641025641026;
	z2.diff(5);

	Dual x3 = 0.568064746731156;
	x3.diff(6);
	Dual y3 = 0.822983865893656;
	y3.diff(7);
	Dual z3 = 0.051282051282051;
	z3.diff(8);

	auto k = compute_curvature(x1, y1, z1, x2, y2, z2, x3, y3, z3);

	double dk_dx1 = k.d(0);
	auto k_x1 = [&](double x) {
		return compute_curvature(x, y1.x(), z1.x(), x2.x(), y2.x(), z2.x(), x3.x(), y3.x(), z3.x());
	};
	double dk_dx1_numeric = (k_x1(x1.x() + h) - k_x1(x1.x() - h)) / (2 * h);
	CAPTURE(dk_dx1);
	CAPTURE(dk_dx1_numeric);
	CHECK(std::abs(dk_dx1 - dk_dx1_numeric) <= 1e-4);

	double dk_dy3 = k.d(7);
	auto k_y3 = [&](double y) {
		return compute_curvature(x1.x(), y1.x(), z1.x(), x2.x(), y2.x(), z2.x(), x3.x(), y, z3.x());
	};
	double dk_dy3_numeric = (k_y3(y3.x() + h) - k_y3(y3.x() - h)) / (2 * h);
	CAPTURE(dk_dy3);
	CAPTURE(dk_dy3_numeric);
	CHECK(std::abs(dk_dy3 - dk_dy3_numeric) <= 1e-4);
}

TEST_CASE("differentiate_torsion") {
	typedef fadbad::F<double, 12> Dual;
	double h = 1e-3;

	Dual x1 = 1;
	x1.diff(0);
	Dual y1 = 0;
	y1.diff(1);
	Dual z1 = 0;
	z1.diff(2);

	Dual x2 = 0.885456025653210;
	x2.diff(3);
	Dual y2 = 0.464723172043769;
	y2.diff(4);
	Dual z2 = 0.025641025641026;
	z2.diff(5);

	Dual x3 = 0.568064746731156;
	x3.diff(6);
	Dual y3 = 0.822983865893656;
	y3.diff(7);
	Dual z3 = 0.051282051282051;
	z3.diff(8);

	Dual x4 = 0.120536680255323;
	x4.diff(9);
	Dual y4 = 0.992708874098054;
	y4.diff(10);
	Dual z4 = 0.076923076923077;
	z4.diff(11);

	auto t = compute_torsion(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

	double dt_dx1 = t.d(0);
	auto t_x1 = [&](double x) {
		return compute_torsion(x,
		                       y1.x(),
		                       z1.x(),
		                       x2.x(),
		                       y2.x(),
		                       z2.x(),
		                       x3.x(),
		                       y3.x(),
		                       z3.x(),
		                       x4.x(),
		                       y4.x(),
		                       z4.x());
	};
	double dt_dx1_numeric = (t_x1(x1.x() + h) - t_x1(x1.x() - h)) / (2 * h);
	CAPTURE(dt_dx1);
	CAPTURE(dt_dx1_numeric);
	CHECK(std::abs(dt_dx1 - dt_dx1_numeric) <= 1e-4);

	double dt_dy4 = t.d(10);
	auto t_y4 = [&](double y) {
		return compute_torsion(x1.x(),
		                       y1.x(),
		                       z1.x(),
		                       x2.x(),
		                       y2.x(),
		                       z2.x(),
		                       x3.x(),
		                       y3.x(),
		                       z3.x(),
		                       x4.x(),
		                       y,
		                       z4.x());
	};
	double dt_dy4_numeric = (t_y4(y4.x() + h) - t_y4(y4.x() - h)) / (2 * h);
	CAPTURE(dt_dy4);
	CAPTURE(dt_dy4_numeric);
	CHECK(std::abs(dt_dy4 - dt_dy4_numeric) <= 1e-4);
}
