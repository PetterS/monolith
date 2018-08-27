#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
using namespace std;

#include <catch.hpp>

#include <minimum/constrained/successive.h>
#include <minimum/core/numeric.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/solver.h>
using minimum::core::relative_error;
using minimum::nonlinear::AutoDiffTerm;
using minimum::nonlinear::ConstrainedFunction;

minimum::constrained::SuccessiveLinearProgrammingSolver::Result solve(
    ConstrainedFunction* function) {
	minimum::constrained::SuccessiveLinearProgrammingSolver solver;
	solver.log_function = [](string_view s) { cerr << s << endl; };
	return solver.solve(function);
}

//
// Ex14_1_3
// http://www.mat.univie.ac.at/~neum/glopt/coconut/Benchmark/Library1/ex14_1_3.mod
//
struct Ex14_1_3_Objective {
	template <typename R>
	R operator()(const R* const x) const {
		return x[2];
	}
};

struct Ex14_1_3_Constraint_e2 {
	template <typename R>
	R operator()(const R* const x) const {
		return 10000 * x[0] * x[1] - x[2] - 1.0;
	}
};

struct Ex14_1_3_Constraint_e3 {
	template <typename R>
	R operator()(const R* const x) const {
		return -10000 * x[0] * x[1] - x[2] + 1.0;
	}
};

struct Ex14_1_3_Constraint_e4 {
	template <typename R>
	R operator()(const R* const x) const {
		return exp(-x[0]) + exp(-x[1]) - x[2] - 1.001;
	}
};

struct Ex14_1_3_Constraint_e5 {
	template <typename R>
	R operator()(const R* const x) const {
		return -exp(-x[0]) - exp(-x[1]) - x[2] + 1.001;
	}
};

struct Ex14_1_3_Constraint_x0_lb {
	template <typename R>
	R operator()(const R* const x) const {
		// x[0] >= 5.49E-6.
		return -x[0] + 5.49E-6;
	}
};

struct Ex14_1_3_Constraint_x0_ub {
	template <typename R>
	R operator()(const R* const x) const {
		return x[0] - 4.553;
	}
};

struct Ex14_1_3_Constraint_x1_lb {
	template <typename R>
	R operator()(const R* const x) const {
		// x[1] >= 0.0021961.
		return -x[1] + 0.0021961;
	}
};

struct Ex14_1_3_Constraint_x1_ub {
	template <typename R>
	R operator()(const R* const x) const {
		return x[1] - 18.21;
	}
};

TEST_CASE("Ex14_1_3") {
	ConstrainedFunction function;
	function.max_number_of_iterations = 1000;

	vector<double> x(3);
	x[0] = 5.49E-6;
	x[1] = 0.0021961;
	x[2] = 0;  // Unspecified by problem.

	function.add_term<AutoDiffTerm<Ex14_1_3_Objective, 3>>(&x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_e2, 3>>("e2", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_e3, 3>>("e3", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_e4, 3>>("e4", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_e5, 3>>("e5", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_x0_lb, 3>>("x0_lb", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_x0_ub, 3>>("x0_ub", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_x1_lb, 3>>("x1_lb", &x[0]);
	function.add_constraint_term<AutoDiffTerm<Ex14_1_3_Constraint_x1_ub, 3>>("x1_ub", &x[0]);

	solve(&function);
	CHECK(x[0] == Approx(0.0000145067));
	CHECK(x[1] == Approx(6.8933528699));
	CHECK(abs(x[2]) <= 1e-9);
	CHECK(function.is_feasible());
	CHECK(abs(function.objective().evaluate()) <= 1e-9);
}

//  (x-a)^2 + (y-b)^2 ≤ 1
class UnitDisc {
   public:
	UnitDisc(double a_, double b_) : a{a_}, b{b_} {}

	template <typename R>
	R operator()(const R* const x) const {
		auto dx = x[0] - a;
		auto dy = x[1] - b;
		return dx * dx + dy * dy - 1;
	}

	template <typename R>
	R operator()(const R* const x, const R* const y) const {
		auto dx = x[0] - a;
		auto dy = y[0] - b;
		return dx * dx + dy * dy - 1;
	}

   private:
	double a, b;
};

// Test in which two discs intersect in a single point.
// The objective function is all zeros.
TEST_CASE("Two_unit_discs") {
	vector<pair<double, double>> start_values = {{3., 4.}, {-10., 10.}, {0., 0.}};

	for (auto start : start_values) {
		ConstrainedFunction function;

		// Infeasible start.
		vector<double> x(3);
		x[0] = start.first;
		x[1] = start.second;
		CAPTURE(x[0]);
		CAPTURE(x[1]);

		function.add_constraint_term(
		    "Disc(0, 0)", make_shared<AutoDiffTerm<UnitDisc, 2>>(0, 0), &x[0]);
		function.add_constraint_term(
		    "Disc(2, 0)", make_shared<AutoDiffTerm<UnitDisc, 2>>(2, 0), &x[0]);

		solve(&function);

		CAPTURE(x[0]);
		CAPTURE(x[1]);
		CHECK(abs(x[0] - 1.0) <= 1e-5);
		CHECK(abs(x[1] - 0.0) <= 1e-5);
		CHECK(function.is_feasible());
	}
}

TEST_CASE("Two_unit_discs_two_variables") {
	vector<pair<double, double>> start_values = {{3., 4.}, {-10., 10.}, {0., 0.}};

	for (auto start : start_values) {
		ConstrainedFunction function;

		// Infeasible start.
		double x = start.first;
		double y = start.second;
		CAPTURE(x);
		CAPTURE(y);

		function.add_constraint_term(
		    "Disc(0, 0)", make_shared<AutoDiffTerm<UnitDisc, 1, 1>>(0, 0), &x, &y);
		function.add_constraint_term(
		    "Disc(2, 0)", make_shared<AutoDiffTerm<UnitDisc, 1, 1>>(2, 0), &x, &y);

		solve(&function);

		CAPTURE(x);
		CAPTURE(y);
		CHECK(abs(x - 1.0) <= 1e-5);
		CHECK(abs(y - 0.0) <= 1e-5);
		CHECK(function.is_feasible());
	}
}

// a·x + b·y ≤ c
class Halfplane {
   public:
	Halfplane(double a_, double b_, double c_) : a{a_}, b{b_}, c{c_} {}

	template <typename R>
	R operator()(const R* const x) const {
		return a * x[0] + b * x[1] - c;
	}

   private:
	double a, b, c;
};

// Test in which the unit circle intersects a line.
// The objective function is all zeros.
TEST_CASE("One_unit_circle_one_line") {
	vector<pair<double, double>> start_values = {{3., 4.}, {10., 20.}, {0., 0.1}};

	for (auto start : start_values) {
		ConstrainedFunction function;

		// Infeasible start.
		vector<double> x(3);
		x[0] = start.first;
		x[1] = start.second;
		CAPTURE(x[0]);
		CAPTURE(x[1]);

		function.add_equality_constraint_term(
		    "Disc(0, 0)", make_shared<AutoDiffTerm<UnitDisc, 2>>(0, 0), &x[0]);
		// The second equality constraint is represented as two inequality constraints
		// just to test that it works.
		function.add_constraint_term(
		    "Plane1", make_shared<AutoDiffTerm<Halfplane, 2>>(1, 1, 0), &x[0]);
		function.add_constraint_term(
		    "Plane2", make_shared<AutoDiffTerm<Halfplane, 2>>(-1, -1, 0), &x[0]);

		solve(&function);

		CAPTURE(x[0]);
		CAPTURE(x[1]);
		CHECK(function.is_feasible());
	}
}
