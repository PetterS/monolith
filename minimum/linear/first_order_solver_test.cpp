#include <cmath>
#include <fstream>

#include <catch.hpp>

#include <minimum/linear/data/util.h>
#include <minimum/linear/first_order_solver.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
using namespace minimum::linear;

TEST_CASE("Strange_convergence") {
	using namespace Eigen;
	using namespace std;

	int m = 4;
	int n = 8;
	MatrixXd Adense(m, n);
	Adense.row(0) << 1, 5, 6, 1, -5, 1, 6, 0;
	Adense.row(1) << 2, 5, 5, -1, 5, 1, 6, 0;
	Adense.row(2) << 3, 5, -4, -1, 5, 2, 0, 9;
	Adense.row(3) << 4, 5, 1, 1, -5, 3, 0, 9;
	minimum::linear::Matrix A;
	A = Adense.sparseView();

	VectorXd x(n);
	x << 1, 1, 0, 0, 1, 1, 0, 0;

	VectorXd b(m);
	b = A * x;
	x.setZero();

	VectorXd c(n);
	c << 1, -1, 2, -3, -5, 7, 9, -6;

	VectorXd y(m);
	y.setZero();

	VectorXd lb(n);
	lb.setZero();

	VectorXd ub(n);
	ub.setConstant(1.0);

	FirstOrderOptions options;
	options.rescale_c = false;

	// From my first implementation, which obtains similar
	// convergence plots as the paper.
	const vector<double> ground_truth = {
	    -5.550000,  -11.359921, -14.174585, -14.017215, -13.926363, -13.702674, -13.234735,
	    -12.690176, -12.339068, -12.289323, -12.469620, -12.722200, -12.900454, -12.898698,
	    -12.640000, -12.260535, -11.814795, -11.328464, -10.804136, -10.233576};

	vector<LinearConstraintType> constraint_types(m, LinearConstraintType::Equality);

	for (size_t i = 0; i < ground_truth.size(); ++i) {
		CAPTURE(i);
		CAPTURE(c.dot(x));
		CAPTURE(ground_truth[i]);

		options.maximum_iterations = i + 1;
		x.setZero();
		y.setZero();
		first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, constraint_types, options);
		REQUIRE(abs(c.dot(x) - ground_truth[i]) < 1e-6);
	}

	const vector<double> ground_truth_admm = {-14.5663, -17.8405, -16.0170, -15.6621, -14.4768,
	                                          -14.0220, -13.8674, -13.8093, -13.2956, -12.4653,
	                                          -11.3204, -9.8875,  -9.2396,  -9.5241,  -9.1223,
	                                          -8.7106,  -8.2899,  -7.8555,  -7.4031,  -6.9324};

	for (size_t i = 0; i < ground_truth_admm.size(); ++i) {
		CAPTURE(i);
		CAPTURE(c.dot(x));
		CAPTURE(ground_truth_admm[i]);

		options.maximum_iterations = i + 1;
		x.setZero();
		first_order_admm_solve(&x, c, lb, ub, A, b, options);
		REQUIRE((abs(c.dot(x) - ground_truth_admm[i]) / abs(c.dot(x))) < 1e-5);
	}

	options.maximum_iterations = 20000;
	options.print_interval = 1000;
	options.log_function = [](const string& str) { cerr << str << endl; };
	x.setZero();
	y.setZero();
	first_order_primal_dual_solve(&x, &y, c, lb, ub, A, b, constraint_types, options);

	cerr << "x     = " << x.transpose() << endl;
	cerr << "<c|x> = " << c.dot(x) << endl;

	CHECK(abs(c.dot(x) - 2.0) <= 1e-6);

	options.maximum_iterations = 2000;
	options.print_interval = 100;
	options.tolerance = 0;
	x.setZero();
	first_order_admm_solve(&x, c, lb, ub, A, b, options);

	cerr << "x     = " << x.transpose() << endl;
	cerr << "<c|x> = " << c.dot(x) << endl;

	CHECK(abs(c.dot(x) - 2.0) <= 1e-6);
}

TEST_CASE("FirstOrderProblem/tiny") {
	using namespace std;

	IP problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.add_bounds(0.0, x, 1.0);
	problem.add_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y == 1.0);
	problem.add_objective(-2.0 * x - 1.0 * y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 1;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&problem)->get());
	CHECK(x.value() >= 0.999);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/tiny/inequality1") {
	using namespace std;

	IP problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.add_bounds(0.0, x, 1.0);
	problem.add_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y <= 1.0);
	problem.add_objective(-2.0 * x - 1.0 * y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 1;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&problem)->get());
	CHECK(x.value() >= 0.999);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/tiny/inequality2") {
	using namespace std;

	IP problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.add_bounds(0.0, x, 1.0);
	problem.add_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y <= 1.0);
	problem.add_objective(2.0 * x + 1.0 * y);

	FirstOrderOptions options;
	options.maximum_iterations = 200;
	options.print_interval = 10;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&problem)->get());
	CHECK(x.value() <= 0.001);
	CHECK(y.value() <= 0.001);
}

TEST_CASE("FirstOrderProblem/no_integer_variables") {
	IP problem;
	problem.add_variable(IP::Integer);
	PrimalDualSolver primal_dual_solver;
	CHECK_THROWS(primal_dual_solver.solutions(&problem)->get());
}

TEST_CASE("FirstOrderProblem/inequalities") {
	using namespace std;

	IP problem;
	auto x1 = problem.add_variable(IP::Real);
	auto x2 = problem.add_variable(IP::Real);
	auto x3 = problem.add_variable(IP::Real);
	auto x4 = problem.add_variable(IP::Real);
	auto x5 = problem.add_variable(IP::Real);

	problem.add_bounds(-10.0, x1, 10.0);
	problem.add_bounds(-10.0, x2, 10.0);
	problem.add_bounds(-10.0, x3, 10.0);
	problem.add_bounds(-10.0, x4, 10.0);
	problem.add_bounds(-10.0, x5, 10.0);

	problem.add_constraint(x1 + 2 * x5 <= 5.0);
	problem.add_constraint(x2 + x4 <= 7.0);
	problem.add_constraint(2 * x3 + x4 <= 4.0);
	problem.add_constraint(x1 + x4 <= 5.0);
	problem.add_constraint(x3 + 2 * x5 <= 4.0);
	problem.add_constraint(x1 + x2 <= 6.0);

	auto obj = -2.0 * x1 - 1.0 * x2 - 3.0 * x3 - 1.0 * x4 - 2.0 * x5;
	problem.add_objective(obj);

	solve(&problem);
	auto optimal_objective = obj.value();
	cerr << "Optimal value: " << optimal_objective << endl;

	FirstOrderOptions options;
	options.maximum_iterations = 20000;
	options.print_interval = 10;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&problem)->get());
	CHECK(abs(obj.value() - optimal_objective) <= 1e-4);
}

TEST_CASE("FirstOrderProblem/inequalities2") {
	using namespace std;

	IP problem;
	auto x1 = problem.add_variable(IP::Real);
	auto x2 = problem.add_variable(IP::Real);
	auto x3 = problem.add_variable(IP::Real);
	auto x4 = problem.add_variable(IP::Real);
	auto x5 = problem.add_variable(IP::Real);

	problem.add_bounds(-10.0, x1, 10.0);
	problem.add_bounds(-10.0, x2, 10.0);
	problem.add_bounds(-10.0, x3, 10.0);
	problem.add_bounds(-10.0, x4, 10.0);
	problem.add_bounds(-10.0, x5, 10.0);

	problem.add_constraint(-x1 - 2 * x5 >= -5.0);
	problem.add_constraint(x2 + x4 <= 7.0);
	problem.add_constraint(2 * x3 + x4 <= 4.0);
	problem.add_constraint(-x1 - x4 >= -5.0);
	problem.add_constraint(-x3 - 2 * x5 >= -4.0);
	problem.add_constraint(x1 + x2 <= 6.0);

	auto obj = -2.0 * x1 - 1.0 * x2 - 3.0 * x3 - 1.0 * x4 - 2.0 * x5;
	problem.add_objective(obj);

	solve(&problem);
	auto optimal_objective = obj.value();
	cerr << "Optimal value: " << optimal_objective << endl;

	FirstOrderOptions options;
	options.maximum_iterations = 20000;
	options.print_interval = 10;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&problem)->get());
	CHECK(abs(obj.value() - optimal_objective) <= 1e-4);
}

TEST_CASE("FirstOrderProblem/infeasible") {
	using namespace std;

	IP problem;
	auto x = problem.add_variable(IP::Real);
	auto y = problem.add_variable(IP::Real);

	problem.add_bounds(0.0, x, 1.0);
	problem.add_bounds(0.0, y, 1.0);
	problem.add_constraint(x + y == 3.0);
	problem.add_objective(-2.0 * x - 1.0 * y);

	FirstOrderOptions options;
	options.maximum_iterations = 2000;
	options.print_interval = 100;
	options.log_function = [](const string& str) { cerr << str << endl; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK_FALSE(primal_dual_solver.solutions(&problem)->get());
}

TEST_CASE("sudoku_warm_start") {
	int n = 3;
	IP ip;
	auto P = ip.add_cube(n * n, n * n, n * n, IP::Real);

	const char* given[] = {
	    "1**"
	    "***"
	    "7**",
	    "**7"
	    "1*9"
	    "***",
	    "68*"
	    "*7*"
	    "***",

	    "**1"
	    "*9*"
	    "6**",
	    "***"
	    "3**"
	    "*2*",
	    "*4*"
	    "***"
	    "**3",

	    "**8"
	    "*6*"
	    "1**",
	    "5**"
	    "***"
	    "*4*",
	    "***"
	    "**2"
	    "**5"};

	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			for (int k = 0; k < n * n; ++k) {
				ip.add_bounds(0, P[i][j][k], 1);
			}
		}
	}

	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			if (given[i][j] != '*') {
				int k = given[i][j] - '1';
				ip.add_constraint(P[i][j][k] == 1);
			}
		}
	}

	// Exactly one indicator equal to 1.
	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			Sum k_sum;
			for (int k = 0; k < n * n; ++k) {
				k_sum += P[i][j][k];
			}
			ip.add_constraint(k_sum == 1);
		}
	}

	// All rows have every number.
	for (int i = 0; i < n * n; ++i) {
		for (int k = 0; k < n * n; ++k) {
			Sum row_k_sum;
			for (int j = 0; j < n * n; ++j) {
				row_k_sum += P[i][j][k];
			}
			ip.add_constraint(row_k_sum == 1);
		}
	}

	// All columns have every number.
	for (int j = 0; j < n * n; ++j) {
		for (int k = 0; k < n * n; ++k) {
			Sum col_k_sum;
			for (int i = 0; i < n * n; ++i) {
				col_k_sum += P[i][j][k];
			}
			ip.add_constraint(col_k_sum == 1);
		}
	}

	// The n*n subsquares have every number.
	for (int i1 = 0; i1 < n; ++i1) {
		for (int j1 = 0; j1 < n; ++j1) {
			for (int k = 0; k < n * n; ++k) {
				Sum square_k_sum;
				for (int i2 = 0; i2 < n; ++i2) {
					for (int j2 = 0; j2 < n; ++j2) {
						square_k_sum += P[n * i1 + i2][n * j1 + j2][k];
					}
				}
				ip.add_constraint(square_k_sum == 1);
			}
		}
	}

	FirstOrderOptions options;
	options.maximum_iterations = 2000;
	options.log_function = [](const std::string& str) { std::cerr << str << "\n"; };
	options.rescale_c = false;

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	CHECK(primal_dual_solver.solutions(&ip)->get());
	CHECK(ip.is_feasible(1e-6));
	CHECK(ip.is_dual_feasible(1e-3));

	// Resolve with warm start.
	primal_dual_solver.options().maximum_iterations = 10;
	CHECK(primal_dual_solver.solutions(&ip)->get());
	CHECK(ip.is_feasible(1e-6));
}

double relative_error(double ground_truth, double calculated) {
	double err = std::abs(ground_truth - calculated);
	return err / std::abs(ground_truth);
}

void netlib_test(const std::string& name, bool rescale_c) {
	std::ifstream fin(data::get_directory() + "/netlib/" + name + ".SIF");
	auto ip = read_MPS(fin);

	GlpkSolver solver;
	solver.set_silent(true);
	auto exact_solution = solver.solutions(ip.get());
	REQUIRE(exact_solution->get());
	double ground_truth = ip->get_entire_objective();
	CHECK(ip->is_feasible());
	CHECK(ip->is_dual_feasible());

	std::stringstream first_order_log;
	FirstOrderOptions options;
	options.maximum_iterations = 1'000'000;
	options.tolerance = 1e-6;
	options.rescale_c = rescale_c;
	options.log_function = [&](const std::string& str) { first_order_log << str << "\n"; };

	PrimalDualSolver primal_dual_solver;
	primal_dual_solver.options() = options;
	bool result = primal_dual_solver.solutions(ip.get())->get();
	CAPTURE(first_order_log.str());
	REQUIRE(result);
	double calulated_first_order_objective = ip->get_entire_objective();
	CHECK(relative_error(ground_truth, calulated_first_order_objective) < 1e-2);
	CHECK(ip->is_feasible(1e-3));
	CHECK(ip->is_dual_feasible(1e-3));
}

TEST_CASE("netlib-AFIRO") { netlib_test("AFIRO", false); }

TEST_CASE("netlib-AFIRO-rescale") { netlib_test("AFIRO", true); }

TEST_CASE("netlib-SC105") { netlib_test("SC105", false); }

TEST_CASE("netlib-SC50A") { netlib_test("SC50A", false); }
