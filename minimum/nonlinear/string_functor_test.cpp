#include <catch.hpp>
using Catch::Detail::Approx;

#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/dynamic_auto_diff_term.h>
#include <minimum/nonlinear/string_functor.h>
using namespace minimum::nonlinear;

TEST_CASE("1_var") {
	AutoDiffTerm<StringFunctor1, 2> term("(x + 2*y)^2", std::vector<std::string>{"x", "y"});
	double vars[] = {5, 1};
	std::vector<double*> all_vars = {vars};

	CHECK(term.evaluate(all_vars.data()) == Approx(49));

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(2));
	CHECK(term.evaluate(all_vars.data(), &gradient) == Approx(49));
	CHECK(gradient[0][0] == Approx(14));
	CHECK(gradient[0][1] == Approx(28));
}

TEST_CASE("1_var_dynamic") {
	DynamicAutoDiffTerm<StringFunctor1> term(
	    {2}, "(x + 2*y)^2", std::vector<std::string>{"x", "y"});
	double vars[] = {5, 1};
	std::vector<double*> all_vars = {vars};

	CHECK(term.evaluate(all_vars.data()) == Approx(49));

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(2));
	CHECK(term.evaluate(all_vars.data(), &gradient) == Approx(49));
	CHECK(gradient[0][0] == Approx(14));
	CHECK(gradient[0][1] == Approx(28));
}

TEST_CASE("n_var_dynamic") {
	DynamicAutoDiffTerm<StringFunctorN> term(
	    {1, 1}, "(x + 2*y)^2", std::vector<std::string>{"x", "y"});
	double x = 5;
	double y = 1;
	std::vector<double*> all_vars = {&x, &y};

	CHECK(term.evaluate(all_vars.data()) == Approx(49));

	std::vector<Eigen::VectorXd> gradient;
	gradient.push_back(Eigen::VectorXd(1));
	gradient.push_back(Eigen::VectorXd(1));
	CHECK(term.evaluate(all_vars.data(), &gradient) == Approx(49));
	CHECK(gradient[0][0] == Approx(14));
	CHECK(gradient[1][0] == Approx(28));
}
