
#include <algorithm>
#include <array>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/first_order_solver.h>
#include <minimum/linear/scs.h>
using namespace minimum::core;
using namespace minimum::linear;

template <typename T, class Container>
vector<T> subvector(const Container& v, size_t n) {
	if (n >= v.size()) {
		return {v.begin(), v.end()};
	} else {
		return {v.begin(), v.begin() + n};
	}
}

int main_program() {
	// constexpr int n = 50'000;
	constexpr int n = 10'000;
	constexpr int ndim = 600;
	using Point = array<double, ndim>;

	mt19937_64 engine(1u);
	uniform_real_distribution<double> coef(-1, 1);
	normal_distribution<double> norm(0, 0.05);

	double b = coef(engine);
	Point a;
	generate(a.begin(), a.end(), bind(coef, engine));

	vector<pair<Point, double>> points(n);

	Timer t("Creating data points");
	for (int i : range(n)) {
		generate(points[i].first.begin(), points[i].first.end(), bind(coef, engine));
		points[i].second =
		    inner_product(points[i].first.begin(), points[i].first.end(), a.begin(), b)
		    + 0.01 * norm(engine);
	}
	t.OK();

	IP ip;
	auto a_var = ip.add_vector(ndim, IP::Real);
	auto b_var = ip.add_variable(IP::Real);

	t.next("Building LP");
	for (int i : range(n)) {
		Sum value = inner_product(
		    points[i].first.begin(), points[i].first.end(), a_var.begin(), Sum(b_var));
		ip.add_objective(ip.abs(value - points[i].second));
	}
	t.OK();

	cerr << ip.get_number_of_variables() << " variables and " << ip.get().constraint_size()
	     << " contraints.\n";

	t.next("Solving with primal-dual");
	ip.clear_solution();
	PrimalDualSolver primal_dual;
	stringstream log;
	primal_dual.options().log_function = [&log](string s) { log << s << endl; };
	primal_dual.options().maximum_iterations = 100'000;
	primal_dual.options().tolerance = 1e-5;
	if (!primal_dual.solutions(&ip)->get()) {
		t.fail();
	} else {
		t.OK();
	}
	double primal_dual_objective = ip.get_entire_objective();
	cerr << log.str();
	cerr << "a  = " << to_string(subvector<double>(a, 20)) << " b  = " << b << "\n";
	cerr << "a* = " << to_string(subvector<Variable>(a_var, 20)) << " b* = " << b_var << "\n";
	cerr << "Objective: " << primal_dual_objective << "\n";

	t.next("Solving with SCS");
	ip.clear_solution();
	ScsSolver scs;
	scs.set_max_iterations(100'000);
	auto scs_solution = scs.solutions(&ip);
	if (!scs_solution->get()) {
		t.fail();
	} else {
		t.OK();
	}
	cerr << scs_solution->log() << "\n";
	double scs_objective = ip.get_entire_objective();
	cerr << "a  = " << to_string(subvector<double>(a, 20)) << " b  = " << b << "\n";
	cerr << "a* = " << to_string(subvector<Variable>(a_var, 20)) << " b* = " << b_var << "\n";
	cerr << "Objective: " << scs_objective << "\n";

	return 0;
}

int main(int argc, char* argv[]) {
	try {
		gflags::ParseCommandLineFlags(&argc, &argv, true);
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
