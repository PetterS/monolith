// Petter Strandmark.

#include <cstddef>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <random>

#include <minimum/core/check.h>
#include <minimum/core/time.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/solver.h>
using minimum::core::check;
using minimum::core::wall_time;

using namespace minimum::nonlinear;

struct LinearObjective {
	double c;
	LinearObjective(double c) { this->c = c; }

	template <typename R>
	R operator()(const R* const x) const {
		return c * x[0];
	}
};

// Barrier for aTx <= b.
class LogBarrierTerm : public Term {
   public:
	std::vector<double> a;
	double b;
	double* mu;
	LogBarrierTerm(const std::vector<double>& a, double b, double* mu) {
		this->a = a;
		this->b = b;
		this->mu = mu;
	}

	virtual int number_of_variables() const override { return int(a.size()); }

	virtual int variable_dimension(int var) const override { return 1; }

	virtual double evaluate(double* const* const x) const override {
		double constraint = -b;
		for (size_t i = 0; i < a.size(); ++i) {
			constraint += a[i] * x[i][0];
		}

		// constraint <= 0 is converted to a barrier
		// - mu * log(-constraint).
		return -(*mu) * std::log(-constraint);
	}

	virtual double evaluate(double* const* const x,
	                        std::vector<Eigen::VectorXd>* gradient) const override {
		double constraint = -b;
		for (size_t i = 0; i < a.size(); ++i) {
			constraint += a[i] * x[i][0];
		}

		auto inside_log = -constraint;
		auto minus_mu_inv_inside_log = -(*mu) * 1.0 / inside_log;
		for (size_t i = 0; i < a.size(); ++i) {
			(*gradient)[i](0) = minus_mu_inv_inside_log * (-a[i]);
		}

		// constraint <= 0 is converted to a barrier
		// - mu * log(-constraint).
		return -(*mu) * std::log(inside_log);
	}

	virtual double evaluate(double* const* const variables,
	                        std::vector<Eigen::VectorXd>* gradient,
	                        std::vector<std::vector<Eigen::MatrixXd>>* hessian) const override {
		check(false, "Hessians not supported.");
		return 0.0;
	}
};

int main_function() {
	std::mt19937 prng(0);
	std::normal_distribution<double> normal;
	auto randn = std::bind(normal, prng);

	const int n = 10;

	// Variables.
	std::vector<double> x(n, 0.0);
	std::vector<double*> all_x_pointers;
	for (auto& single_x : x) {
		all_x_pointers.push_back(&single_x);
	}

	// Barrier parameter.
	double mu = 1.0;

	// Objective function.
	Function f;

	// Generate random objective vector.
	std::vector<double> c(n);
	for (size_t i = 0; i < n; ++i) {
		c[i] = randn();
		f.add_term(std::make_shared<AutoDiffTerm<LinearObjective, 1>>(c[i]), &x[i]);
	}

	// sum x_i <= 10
	double b = 10;
	std::vector<double> a1(n, 1.0);
	f.add_term(std::make_shared<LogBarrierTerm>(a1, b, &mu), all_x_pointers);

	//  sum x_i >= -10  <=>
	// -sum x_i <= 10
	std::vector<double> a2(n, -1.0);
	b = 10;
	f.add_term(std::make_shared<LogBarrierTerm>(a2, b, &mu), all_x_pointers);

	// Add barriers for individual scalars.
	//
	//		-100 <= x[i] <= 100
	//
	std::vector<double> a3(1, 1.0);
	auto individual_barrier_high = std::make_shared<LogBarrierTerm>(a3, 100, &mu);
	std::vector<double> a4(1, -1.0);
	auto individual_barrier_low = std::make_shared<LogBarrierTerm>(a4, 100, &mu);
	for (size_t i = 0; i < n; ++i) {
		f.add_term(individual_barrier_high, &x[i]);
		f.add_term(individual_barrier_low, &x[i]);
	}

	LBFGSSolver solver;
	solver.maximum_iterations = 100;
	solver.log_function = nullptr;
	SolverResults results;

	double total_time = 0.0;
	for (int iter = 1; iter <= 8; ++iter) {
		auto start_time = wall_time();
		solver.solve(f, &results);

		double sumx = 0.0;
		double cTx = 0.0;
		double minx = std::numeric_limits<double>::infinity();
		double maxx = -minx;
		for (size_t i = 0; i < n; ++i) {
			cTx += x[i] * c[i];
			sumx += x[i];
			minx = std::min(x[i], minx);
			maxx = std::max(x[i], maxx);
		}

		total_time += wall_time() - start_time;
		std::printf("mu=%.2e, f=%+.4e, cTx=%+.4e, sum(x)=%+.4e, (min,max)(x)=(%.3e, %3e)\n",
		            mu,
		            f.evaluate(),
		            cTx,
		            sumx,
		            minx,
		            maxx);

		mu /= 10.0;
	}

	// Check that all constraints are satisfied.
	double sum = 0;
	double eps = 1e-8;
	for (size_t i = 0; i < n; ++i) {
		minimum_core_assert(-100 - eps <= x[i] && x[i] <= 100 + eps, "x[", i, "] = ", x[i]);
		sum += x[i];
		if (i < 10) {
			std::cout << "c[" << i << "] = " << c[i] << ", x[" << i << "] = " << x[i] << std::endl;
		}
	}
	minimum_core_assert(-10 - eps <= sum && sum <= 10 + eps, "sum = ", sum);

	std::cout << "Solution to the linear programming problem is cTx = " << f.evaluate() << '\n';
	std::cout << "Elapsed time: " << total_time << " s.\n";

	return 0;
}

int main() {
	try {
		return main_function();
	} catch (std::exception& e) {
		std::cerr << "Exception: " << e.what() << std::endl;
	}
}
