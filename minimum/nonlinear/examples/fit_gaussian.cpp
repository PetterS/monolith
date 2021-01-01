// Petter Strandmark.

#include <functional>
#include <iostream>
#include <random>

#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/interval_term.h>
#include <minimum/nonlinear/solver.h>
#include <minimum/nonlinear/transformations.h>
using minimum::core::from_string;
using minimum::core::wall_time;

using namespace minimum::nonlinear;

// One term in the negative log-likelihood function for
// a one-dimensional Gaussian distribution.
struct NegLogLikelihood {
	double sample;
	NegLogLikelihood(double sample) { this->sample = sample; }

	template <typename R>
	R operator()(const R* const mu, const R* const sigma) const {
		R diff = (*mu - sample) / *sigma;
		return 0.5 * diff * diff + log(*sigma);
	}
};

int main_function(int argc, char* argv[]) {
	using namespace std;

	int n_samples = 1000;
	if (argc > 1) {
		n_samples = from_string<int>(argv[1]);
	}

	mt19937 prng(1u);
	normal_distribution<double> normal;
	auto randn = bind(normal, prng);

	double mu = 5.0;
	double sigma = 3.0;
	cout << "mu = " << mu << ", sigma = " << sigma << endl;

	Function f;
	f.add_variable(&mu, 1);
	f.add_variable_with_change<GreaterThanZero>(&sigma, 1, 1);

	Function f_global;
	f_global.add_variable(&mu, 1);
	f_global.add_variable(&sigma, 1);

	for (int i = 0; i < n_samples; ++i) {
		double sample = sigma * randn() + mu;
		f.add_term(std::make_shared<AutoDiffTerm<NegLogLikelihood, 1, 1>>(sample), &mu, &sigma);
		f_global.add_term(
		    std::make_shared<IntervalTerm<NegLogLikelihood, 1, 1>>(sample), &mu, &sigma);
	}

	mu = 0.0;
	sigma = 1.0;
	LBFGSSolver solver;

	solver.callback_function = [&](const CallbackInformation& information) -> bool {
		f.copy_global_to_user(*information.x);
		cerr << "  -- mu = " << mu << ", sigma = " << sigma << endl;
		return true;
	};

	SolverResults results;
	solver.solve(f, &results);
	cout << "Estimated:" << endl;
	cout << "f = " << f.evaluate() << " mu = " << mu << ", sigma = " << sigma << endl << endl;

	cout << "Perform global optimization? (y/n):" << endl;
	char answer = 'n';
	cin >> answer;

	if (cin && tolower(answer) == 'y') {
		f_global.set_variable_bounds(&mu, {{-10.0, 10.0}});
		f_global.set_variable_bounds(&sigma, {{0.1, 10.0}});

		GlobalSolver global_solver;
		global_solver.maximum_iterations = 10000;
		auto start_time = wall_time();
		auto interval = global_solver.solve_global(f_global, &results);
		auto elapsed_time = wall_time() - start_time;

		cout << "Optimal parameter intervals:" << endl;
		cout << "   mu    = " << interval.at(0) << endl;
		cout << "   sigma = " << interval.at(1) << endl;
		cout << "Elapsed time was " << elapsed_time << " seconds." << endl;
	}
	return 0;
}

int main(int argc, char* argv[]) {
	try {
		return main_function(argc, argv);
	} catch (std::exception& e) {
		std::cerr << "Error: " << e.what() << '\n';
		return 1;
	}
}
