#include <cmath>
#include <limits>
#include <random>

#include <catch.hpp>

#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/solver.h>

using namespace minimum::nonlinear;

TEST_CASE("Rosenbrock") {
	Function f;
	double x[2] = {2, 3};
	f.add_variable(x, 2);
	f.add_term(make_differentiable<2>([](auto x) {
		           auto d0 = x[1] - x[0] * x[0];
		           auto d1 = 1.0 - x[0];
		           return 100 * d0 * d0 + d1 * d1;
	           }),
	           x);

	LBFGSBSolver solver;
	SolverResults results;
	solver.solve(f, &results);
	CHECK(results.exit_success());

	solver.callback_function = [](const CallbackInformation& information) -> bool { return false; };
	x[0] = 2;
	x[1] = 3;
	solver.solve(f, &results);
	CHECK(results.exit_condition == SolverResults::USER_ABORT);
}

TEST_CASE("10_iterations") {
	Function f;
	double x[2] = {2, 3};
	f.add_variable(x, 2);
	f.add_term(make_differentiable<2>([](auto x) {
		           auto d0 = x[1] - x[0] * x[0];
		           auto d1 = 1.0 - x[0];
		           return 100 * d0 * d0 + d1 * d1;
	           }),
	           x);

	LBFGSBSolver solver;
	SolverResults results;
	solver.maximum_iterations = 10;
	solver.solve(f, &results);
	CHECK(results.exit_condition == SolverResults::NO_CONVERGENCE);
}

TEST_CASE("least_squares_01") {
	using namespace std;
	int m = 10;
	int n = 10;

	mt19937_64 rng;
	uniform_real_distribution<double> rand(-5, 5);
	Eigen::MatrixXd A(m, n);
	Eigen::VectorXd b(m);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			A(i, j) = rand(rng);
		}
		b(i) = rand(rng);
	}

	Function f;
	vector<double> x(n, 0);
	f.add_variable(x.data(), n);

	for (int i = 0; i < m; ++i) {
		f.add_term(make_differentiable<10>([&A, &b, i, n](auto x) {
			           typename std::remove_reference<decltype(x[0])>::type sum = 0;
			           for (int j = 0; j < n; ++j) {
				           sum += x[j] * A(i, j);
			           }
			           return (sum - b[i]) * (sum - b[i]);
		           }),
		           x.data());
	}

	f.set_variable_bounds(x.data(), vector<Interval<double>>(n, {0.0, 1.0}));

	LBFGSBSolver solver;
	SolverResults results;
	solver.solve(f, &results);
	CHECK(results.exit_success());

	for (int j = 0; j < n; ++j) {
		cerr << "x[" << j << "] = " << x[j] << "\n";
		CHECK((0 <= x[j] && x[j] <= 1));
	}

	f.set_variable_bounds(x.data(),
	                      vector<Interval<double>>(n, {-1.0, Interval<double>::infinity}));
	solver.solve(f, &results);
	CHECK(results.exit_success());
	for (int j = 0; j < n; ++j) {
		cerr << "x[" << j << "] = " << x[j] << "\n";
		CHECK(-1 <= x[j]);
	}

	f.set_variable_bounds(x.data(),
	                      vector<Interval<double>>(n, {-Interval<double>::infinity, 2.0}));
	solver.solve(f, &results);
	CHECK(results.exit_success());
	for (int j = 0; j < n; ++j) {
		cerr << "x[" << j << "] = " << x[j] << "\n";
		CHECK(x[j] <= 2);
	}
}
