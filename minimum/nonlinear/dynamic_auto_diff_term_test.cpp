#include <random>
#include <sstream>

#include <catch.hpp>

#include <minimum/nonlinear/dynamic_auto_diff_term.h>
#include <minimum/nonlinear/function.h>

using namespace fadbad;
using namespace minimum::nonlinear;

class DistanceFunctor {
   public:
	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const vars) const {
		using std::sqrt;

		R d = 0;
		for (int i = 1; i < dimensions.size(); ++i) {
			R dx = vars[0][0] - vars[i][0];
			R dy = vars[0][1] - vars[i][1];
			d += sqrt(dx * dx + dy * dy);
		}
		return d;
	}
};

TEST_CASE("DynamicAutoDiffTerm/LargeTerms") {
	constexpr int num_variables = 10000;  // Number of 2D variables.
	constexpr int num_terms = 100;
	constexpr double p = 0.01;  // Probability to include a variable.

	std::mt19937_64 engine;
	std::uniform_real_distribution<double> rand(0.0, 1.0);

	Function f;
	std::vector<std::vector<double>> x(num_variables, {0.0, 0.0});
	for (int i = 0; i < num_variables; ++i) {
		f.add_variable(x[i].data(), 2);
		x[i][0] = rand(engine);
		x[i][1] = rand(engine);
	}

	for (int j = 0; j < num_terms; ++j) {
		std::vector<double*> variables_in_term;
		std::vector<int> dimensions;
		for (int i = 0; i < num_variables; ++i) {
			if (rand(engine) < p) {
				// Add this variable to the term.
				variables_in_term.push_back(x[i].data());
				dimensions.push_back(2);
			}
		}
		f.add_term(std::make_shared<DynamicAutoDiffTerm<DistanceFunctor>>(dimensions),
		           variables_in_term);
	}

	Eigen::VectorXd x_vector;
	Eigen::VectorXd gradient(f.get_number_of_scalars());
	f.copy_user_to_global(&x_vector);
	double f_val = f.evaluate(x_vector, &gradient);
	CHECK(f_val == f_val);  // Check for NaN.
}

class MyFunctor1 {
   public:
	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const x) const {
		return sin(x[0][0]) + cos(x[0][1]) + R(1.4) * x[0][0] * x[0][1] + R(1.0);
	}
};

TEST_CASE("DynamicAutoDiffTerm/MyFunctor1", "") {
	DynamicAutoDiffTerm<MyFunctor1> term({2});

	double x[2] = {1.0, 3.0};
	std::vector<double*> variables;
	variables.push_back(x);

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(2));

	std::vector<std::vector<Eigen::MatrixXd>> hessian(1);
	hessian[0].resize(1);
	hessian[0][0].resize(2, 2);

	double value1 = term.evaluate(&variables[0]);
	double value2 = term.evaluate(&variables[0], &gradient);

	// The values must agree.
	CHECK(Approx(value1) == value2);

	// Test function value
	CHECK(Approx(value1) == sin(x[0]) + cos(x[1]) + 1.4 * x[0] * x[1] + 1.0);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == cos(x[0]) + 1.4 * x[1]);
	CHECK(Approx(gradient[0](1)) == -sin(x[1]) + 1.4 * x[0]);

	double value3 = term.evaluate(&variables[0], &gradient, &hessian);
	CHECK(Approx(value1) == value3);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == cos(x[0]) + 1.4 * x[1]);
	CHECK(Approx(gradient[0](1)) == -sin(x[1]) + 1.4 * x[0]);

	// Test Hessian
	CHECK(Approx(hessian[0][0](0, 0)) == -sin(x[0]));
	CHECK(Approx(hessian[0][0](1, 1)) == -cos(x[1]));
	CHECK(Approx(hessian[0][0](0, 1)) == 1.4);
	CHECK(Approx(hessian[0][0](1, 0)) == 1.4);

	CHECK(term.number_of_variables() == 1);
	CHECK(term.variable_dimension(0) == 2);
}

class MyFunctor2 {
   public:
	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const vars) const {
		auto x = vars[0];
		auto y = vars[1];
		return sin(x[0]) + cos(y[0]) + R(1.4) * x[0] * y[0] + R(1.0);
	}
};

TEST_CASE("DynamicAutoDiffTerm/MyFunctor2", "") {
	DynamicAutoDiffTerm<MyFunctor2> term({1, 1});

	double x = 5.3;
	double y = 7.1;
	std::vector<double*> variables;
	variables.push_back(&x);
	variables.push_back(&y);

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(1));
	gradient.push_back(Eigen::VectorXd(1));

	double value = term.evaluate(&variables[0], &gradient);
	double value2 = term.evaluate(&variables[0]);

	// The two values must agree.
	CHECK(Approx(value) == value2);

	// Test function value
	CHECK(Approx(value) == sin(x) + cos(y) + 1.4 * x * y + 1.0);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == cos(x) + 1.4 * y);
	CHECK(Approx(gradient[1](0)) == -sin(y) + 1.4 * x);

	std::vector<std::vector<Eigen::MatrixXd>> hessian(2);
	hessian[0].resize(2);
	hessian[1].resize(2);
	hessian[0][0].resize(1, 1);
	hessian[0][1].resize(1, 1);
	hessian[1][0].resize(1, 1);
	hessian[1][1].resize(1, 1);

	for (auto& g : gradient) {
		g.setZero();
	}
	double value3 = term.evaluate(&variables[0], &gradient, &hessian);
	CHECK(Approx(value3) == value2);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == cos(x) + 1.4 * y);
	CHECK(Approx(gradient[1](0)) == -sin(y) + 1.4 * x);

	// Test Hessian
	CHECK(Approx(hessian[0][0](0, 0)) == -sin(x));
	CHECK(Approx(hessian[1][1](0, 0)) == -cos(y));
	CHECK(Approx(hessian[1][0](0, 0)) == 1.4);
	CHECK(Approx(hessian[0][1](0, 0)) == 1.4);

	CHECK(term.number_of_variables() == 2);
	CHECK(term.variable_dimension(0) == 1);
	CHECK(term.variable_dimension(1) == 1);
}

class MyFunctor3 {
   public:
	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const vars) const {
		auto x = vars[0];
		auto y = vars[1];
		auto z = vars[2];
		return 2.0 * x[0] + 2.0 * y[0] + 3.0 * y[1] + 2.0 * z[0] * z[0] + 3.0 * z[1] * z[1]
		       + 4.0 * z[2] * z[2];
	}
};

TEST_CASE("DynamicAutoDiffTerm/MyFunctor3") {
	DynamicAutoDiffTerm<MyFunctor3> term({1, 2, 3});

	double x[1] = {5.3};
	double y[2] = {7.1, 5.1};
	double z[3] = {9.5, 1.1, 5.2};
	std::vector<double*> variables;
	variables.push_back(x);
	variables.push_back(y);
	variables.push_back(z);

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(1));
	gradient.push_back(Eigen::VectorXd(2));
	gradient.push_back(Eigen::VectorXd(3));

	std::vector<std::vector<Eigen::MatrixXd>> hessian(3);
	hessian[0].resize(3);
	hessian[1].resize(3);
	hessian[2].resize(3);
	hessian[0][0].resize(1, 1);
	hessian[0][1].resize(1, 2);
	hessian[0][2].resize(1, 3);
	hessian[1][0].resize(2, 1);
	hessian[2][0].resize(3, 1);
	hessian[1][1].resize(2, 2);
	hessian[1][2].resize(2, 3);
	hessian[2][1].resize(3, 2);
	hessian[2][2].resize(3, 3);

	double value = term.evaluate(&variables[0], &gradient);
	double value2 = term.evaluate(&variables[0]);

	// The two values must agree.
	CHECK(Approx(value) == value2);

	// Test function value
	CHECK(Approx(value)
	      == (2.0 * x[0] + 2.0 * y[0] + 3.0 * y[1] + 2.0 * z[0] * z[0] + 3.0 * z[1] * z[1]
	          + 4.0 * z[2] * z[2]));

	// Test gradient
	CHECK(Approx(gradient[0](0)) == 2.0);
	CHECK(Approx(gradient[1](0)) == 2.0);
	CHECK(Approx(gradient[1](1)) == 3.0);
	CHECK(Approx(gradient[2](0)) == 2.0 * 2.0 * z[0]);
	CHECK(Approx(gradient[2](1)) == 2.0 * 3.0 * z[1]);
	CHECK(Approx(gradient[2](2)) == 2.0 * 4.0 * z[2]);

	for (auto& g : gradient) {
		g.setZero();
	}
	double value3 = term.evaluate(&variables[0], &gradient, &hessian);
	CHECK(Approx(value3) == value2);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == 2.0);
	CHECK(Approx(gradient[1](0)) == 2.0);
	CHECK(Approx(gradient[1](1)) == 3.0);
	CHECK(Approx(gradient[2](0)) == 2.0 * 2.0 * z[0]);
	CHECK(Approx(gradient[2](1)) == 2.0 * 3.0 * z[1]);
	CHECK(Approx(gradient[2](2)) == 2.0 * 4.0 * z[2]);

	// Test Hessian
	CHECK(Approx(hessian[0][0](0, 0)) == 0.0);

	CHECK(Approx(hessian[1][1](0, 0)) == 0.0);
	CHECK(Approx(hessian[1][1](0, 1)) == 0.0);
	CHECK(Approx(hessian[1][1](1, 0)) == 0.0);
	CHECK(Approx(hessian[1][1](1, 1)) == 0.0);

	CHECK(Approx(hessian[2][2](0, 0)) == 2.0 * 2.0);
	CHECK(Approx(hessian[2][2](1, 1)) == 2.0 * 3.0);
	CHECK(Approx(hessian[2][2](2, 2)) == 2.0 * 4.0);

	CHECK(term.number_of_variables() == 3);
	CHECK(term.variable_dimension(0) == 1);
	CHECK(term.variable_dimension(1) == 2);
	CHECK(term.variable_dimension(2) == 3);
}

class MyFunctor4 {
   public:
	template <typename R>
	R operator()(const std::vector<int>& dimensions, const R* const* const vars) const {
		auto x = vars[0];
		auto y = vars[1];
		auto z = vars[2];
		auto w = vars[3];
		return 2.0 * x[0] * x[0] + 2.0 * y[0] * y[0] + 3.0 * y[1] * y[1] + 2.0 * z[0] * z[0]
		       + 3.0 * z[1] * z[1] + 4.0 * z[2] * z[2] + 2.0 * w[0] * z[0] + 3.0 * w[1] * z[1]
		       + 4.0 * w[2] * z[2] + 5.0 * w[3] * w[3];
	}
};

TEST_CASE("DynamicAutoDiffTerm/MyFunctor4") {
	DynamicAutoDiffTerm<MyFunctor4> term({1, 2, 3, 4});

	double x[1] = {5.3};
	double y[2] = {7.1, 5.1};
	double z[3] = {9.5, 1.1, 5.2};
	double w[4] = {2.1, 7.87, 2.0, -1.9};
	std::vector<double*> variables;
	variables.push_back(x);
	variables.push_back(y);
	variables.push_back(z);
	variables.push_back(w);

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(1));
	gradient.push_back(Eigen::VectorXd(2));
	gradient.push_back(Eigen::VectorXd(3));
	gradient.push_back(Eigen::VectorXd(4));

	std::vector<std::vector<Eigen::MatrixXd>> hessian(4);
	hessian[0].resize(4);
	hessian[1].resize(4);
	hessian[2].resize(4);
	hessian[3].resize(4);

	hessian[0][0].resize(1, 1);
	hessian[0][1].resize(1, 2);
	hessian[0][2].resize(1, 3);
	hessian[0][3].resize(1, 4);
	hessian[1][0].resize(2, 1);
	hessian[2][0].resize(3, 1);
	hessian[3][0].resize(4, 1);

	hessian[1][1].resize(2, 2);
	hessian[1][2].resize(2, 3);
	hessian[1][3].resize(2, 4);
	hessian[2][1].resize(3, 2);
	hessian[3][1].resize(4, 2);

	hessian[2][2].resize(3, 3);
	hessian[2][3].resize(3, 4);
	hessian[3][2].resize(4, 3);

	hessian[3][3].resize(4, 4);

	double value = term.evaluate(&variables[0], &gradient);
	double value2 = term.evaluate(&variables[0]);

	// The two values must agree.
	CHECK(Approx(value) == value2);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == 2.0 * 2.0 * x[0]);

	CHECK(Approx(gradient[1](0)) == 2.0 * 2.0 * y[0]);
	CHECK(Approx(gradient[1](1)) == 2.0 * 3.0 * y[1]);

	CHECK(Approx(gradient[2](0)) == 2.0 * 2.0 * z[0] + 2.0 * w[0]);
	CHECK(Approx(gradient[2](1)) == 2.0 * 3.0 * z[1] + 3.0 * w[1]);
	CHECK(Approx(gradient[2](2)) == 2.0 * 4.0 * z[2] + 4.0 * w[2]);

	CHECK(Approx(gradient[3](0)) == 2.0 * z[0]);
	CHECK(Approx(gradient[3](1)) == 3.0 * z[1]);
	CHECK(Approx(gradient[3](2)) == 4.0 * z[2]);
	CHECK(Approx(gradient[3](3)) == 2.0 * 5.0 * w[3]);

	for (auto& g : gradient) {
		g.setZero();
	}
	double value3 = term.evaluate(&variables[0], &gradient, &hessian);
	CHECK(Approx(value3) == value2);

	// Test gradient
	CHECK(Approx(gradient[0](0)) == 2.0 * 2.0 * x[0]);

	CHECK(Approx(gradient[1](0)) == 2.0 * 2.0 * y[0]);
	CHECK(Approx(gradient[1](1)) == 2.0 * 3.0 * y[1]);

	CHECK(Approx(gradient[2](0)) == 2.0 * 2.0 * z[0] + 2.0 * w[0]);
	CHECK(Approx(gradient[2](1)) == 2.0 * 3.0 * z[1] + 3.0 * w[1]);
	CHECK(Approx(gradient[2](2)) == 2.0 * 4.0 * z[2] + 4.0 * w[2]);

	CHECK(Approx(gradient[3](0)) == 2.0 * z[0]);
	CHECK(Approx(gradient[3](1)) == 3.0 * z[1]);
	CHECK(Approx(gradient[3](2)) == 4.0 * z[2]);
	CHECK(Approx(gradient[3](3)) == 2.0 * 5.0 * w[3]);

	// Test Hessian
	CHECK(Approx(hessian[0][0](0, 0)) == 2.0 * 2.0);

	CHECK(Approx(hessian[1][1](0, 0)) == 2.0 * 2.0);
	CHECK(Approx(hessian[1][1](0, 1)) == 0.0);
	CHECK(Approx(hessian[1][1](1, 0)) == 0.0);
	CHECK(Approx(hessian[1][1](1, 1)) == 2.0 * 3.0);

	CHECK(Approx(hessian[3][3](0, 0)) == 0.0);
	CHECK(Approx(hessian[3][3](1, 1)) == 0.0);
	CHECK(Approx(hessian[3][3](2, 2)) == 0.0);
	CHECK(Approx(hessian[3][3](3, 3)) == 2.0 * 5.0);
	CHECK(Approx(hessian[3][3](0, 1)) == 0.0);
	CHECK(Approx(hessian[3][3](0, 2)) == 0.0);
	CHECK(Approx(hessian[3][3](0, 3)) == 0.0);
	CHECK(Approx(hessian[3][3](1, 0)) == 0.0);
	CHECK(Approx(hessian[3][3](1, 2)) == 0.0);
	CHECK(Approx(hessian[3][3](1, 3)) == 0.0);
	CHECK(Approx(hessian[3][3](3, 2)) == 0.0);

	CHECK(Approx(hessian[2][3](0, 0)) == 2.0);
	CHECK(Approx(hessian[2][3](1, 1)) == 3.0);
	CHECK(Approx(hessian[2][3](2, 2)) == 4.0);
	CHECK(Approx(hessian[2][3](1, 0)) == 0.0);
	CHECK(Approx(hessian[2][3](1, 2)) == 0.0);

	CHECK(Approx(hessian[3][2](0, 0)) == 2.0);
	CHECK(Approx(hessian[3][2](1, 1)) == 3.0);
	CHECK(Approx(hessian[3][2](2, 2)) == 4.0);
	CHECK(Approx(hessian[3][2](1, 0)) == 0.0);
	CHECK(Approx(hessian[3][2](1, 2)) == 0.0);

	CHECK(term.number_of_variables() == 4);
	CHECK(term.variable_dimension(0) == 1);
	CHECK(term.variable_dimension(1) == 2);
	CHECK(term.variable_dimension(2) == 3);
	CHECK(term.variable_dimension(3) == 4);
}
