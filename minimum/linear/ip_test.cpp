#include <type_traits>

#include <catch.hpp>

#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
#include <minimum/linear/testing.h>

using namespace std;
using namespace minimum::linear;

TEST_CASE("empty") {
	IP ip;
	CHECK(solve(&ip));
}

TEST_CASE("solve_boolean_infeasible") {
	IP ip;
	auto x = ip.add_variable();
	auto y = ip.add_variable();
	ip.add_constraint(x + y >= 3);
	CHECK(!solve(&ip));
}

TEST_CASE("solve_boolean_feasible") {
	IP ip;
	auto x = ip.add_variable();
	auto y = ip.add_variable();
	ip.add_constraint(x + y >= 2);
	CHECK(solve(&ip));
}

TEST_CASE("solve_integer_unbounded") {
	IP ip;
	auto x = ip.add_variable(IP::Integer, -1.0);
	auto y = ip.add_variable(IP::Integer, -1.0);
	ip.add_constraint(x + y >= 3);
	CHECK(!solve(&ip));
}

TEST_CASE("solve_integer_infeasible_real_feasible") {
	IP ip;
	auto x = ip.add_variable(IP::Integer, -1.0);
	auto y = ip.add_variable(IP::Integer, -1.0);
	ip.add_constraint(x + y == 1);
	ip.add_constraint(x >= 0.2);
	ip.add_constraint(y >= 0.2);
	CHECK(!solve(&ip));
}

TEST_CASE("solve_real_unbounded") {
	IP ip;
	auto x = ip.add_variable(IP::Real, -1.0);
	auto y = ip.add_variable(IP::Real, -1.0);
	ip.add_constraint(x + y >= 3);
	CHECK(!solve(&ip));
}

TEST_CASE("solve_real_feasible") {
	IP ip;
	auto x = ip.add_variable(IP::Real, -1.0);
	auto y = ip.add_variable(IP::Real, -1.0);
	ip.add_bounds(-100, x, 100);
	ip.add_bounds(-100, y, 100);
	ip.add_constraint(x + y <= 3);
	CHECK(solve(&ip));
}

TEST_CASE("solve_real_infeasible") {
	IP ip;
	auto x = ip.add_variable(IP::Real, -1.0);
	auto y = ip.add_variable(IP::Real, -1.0);
	ip.add_bounds(-100, x, 100);
	ip.add_bounds(-100, y, 100);
	ip.add_constraint(x + y >= 1000);
	CHECK(!solve(&ip));
}

TEST_CASE("guaranteed_infeasible") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	auto y = ip.add_variable(IP::Real);
	ip.add_constraint(x >= 5);
	ip.add_constraint(x <= 4);
	ip.add_constraint(x + y >= 1000);
	CHECK(!solve(&ip));
}

TEST_CASE("guaranteed_infeasible_constraints") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	auto y = ip.add_variable(IP::Real);
	ip.add_bounds(-100, x, 100);
	ip.add_bounds(-100, y, 100);
	ip.add_constraint(5000, x + y, 1000);
	CHECK(!solve(&ip));
}

TEST_CASE("different_solvers") {
	IP ip1, ip2;
	auto x = ip1.add_boolean();
	auto y = ip2.add_boolean();
	CHECK_THROWS(x + y);
	CHECK_THROWS(x - y);
	CHECK_THROWS(x || y);
	CHECK_THROWS(ip1.add_objective(y));
	CHECK_THROWS(ip2.add_bounds(0, x, 2));
	CHECK_THROWS(x == y);
	CHECK_THROWS(x >= y);
	CHECK_THROWS(x <= y);
	CHECK_THROWS(ip1.add_constraint(2.0 * y >= 8));
	CHECK_THROWS(ip2.add_constraint(x + x >= 8));
	solve(&ip1);
	solve(&ip2);
	CHECK_THROWS(ip1.get_solution(y));
	CHECK_THROWS(ip1.get_solution(!y));
}

TEST_CASE("no_copy/move") {
	static_assert(!std::is_move_constructible<IP>::value, "IP should not be movable.");
	static_assert(!std::is_copy_constructible<IP>::value, "IP should not be movable.");
}

TEST_CASE("constraints_compilation") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x <= 0);
	ip.add_constraint(0 <= x);
	ip.add_constraint(x <= y);
	ip.add_constraint(x + y <= 0);
	ip.add_constraint(x + y <= y);
	ip.add_constraint(0 <= x + y);
	ip.add_constraint(y <= x + y);
	ip.add_constraint(x + y <= x + y);

	ip.add_constraint(x >= 0);
	ip.add_constraint(0 >= x);
	ip.add_constraint(x >= y);
	ip.add_constraint(x + y >= 0);
	ip.add_constraint(x + y >= y);
	ip.add_constraint(0 >= x + y);
	ip.add_constraint(y >= x + y);
	ip.add_constraint(x + y >= x + y);

	ip.add_constraint(x == 0);
	ip.add_constraint(0 == x);
	ip.add_constraint(x == y);
	ip.add_constraint(x + y == 0);
	ip.add_constraint(x + y == y);
	ip.add_constraint(0 == x + y);
	ip.add_constraint(y == x + y);
	ip.add_constraint(x + y == x + y);
	SUCCEED();
}

TEST_CASE("constraints_overwritten") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	ip.add_constraint(x <= 100);
	ip.add_constraint(x >= -100);
	ip.add_bounds(-200, x, 200);
	ip.add_objective(x);
	solve(&ip);

	CHECK(Approx(x.value()) == -100.0);
}

TEST_CASE("self_constraint") {
	IP ip;
	auto x = ip.add_variable();
	auto y = ip.add_variable();
	ip.add_constraint(x + y >= 2);
	ip.add_constraint(x <= x);
	ip.add_constraint(y >= y);
	ip.add_constraint(x == x);
	CHECK(solve(&ip));
}

TEST_CASE("grid") {
	IP ip;
	auto x = ip.add_grid(5, 6);
	REQUIRE(x.size() == 5);
	for (auto& row : x) {
		REQUIRE(row.size() == 6);
	}
}

TEST_CASE("constraint_list") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(x == 1 && y == 0 && z == 1);
	REQUIRE(solve(&ip));
	CHECK(x.bool_value());
	CHECK_FALSE(y.bool_value());
	CHECK(z.bool_value());
}

TEST_CASE("bool_constraint") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x);
	ip.add_constraint(y);
	REQUIRE(solve(&ip));
	CHECK(x.bool_value());
	CHECK(y.bool_value());
}

TEST_CASE("basic_ip") {
	IPSolver solver;
	solver.set_silent(true);
	CHECK_NOTHROW(simple_integer_programming_tests(solver));
}

TEST_CASE("basic_lp") {
	IPSolver solver;
	solver.set_silent(true);
	CHECK_NOTHROW(simple_linear_programming_tests(solver));
}

TEST_CASE("add_variable_as_booleans") {
	IP ip;
	auto x = ip.add_variable_as_booleans(-3, 3);
	auto y = ip.add_variable_as_booleans({-3, -2, 3, 2, 1});
	ip.add_objective(x - y);
	CHECK(solve(&ip));

	CHECK(x.value() == -3);
	CHECK(y.value() == 3);
}

TEST_CASE("add_max_consecutive_constraints") {
	IP ip;
	vector<Sum> x;
	for (int i = 1; i <= 10; ++i) {
		x.emplace_back(ip.add_boolean(-1.0));
	}
	ip.add_max_consecutive_constraints(5, x);
	int num_solutions = 0;
	IPSolver solver;
	solver.set_silent(true);
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	REQUIRE(num_solutions == 2);
	CHECK(x[0].value() == 1);
	CHECK(x[1].value() == 1);
	CHECK(x[2].value() == 1);
	CHECK(x[3].value() == 1);
	// CHECK(x[4].value() == ?);
	// CHECK(x[5].value() == ?);
	CHECK(x[6].value() == 1);
	CHECK(x[7].value() == 1);
	CHECK(x[8].value() == 1);
	CHECK(x[9].value() == 1);
}

TEST_CASE("add_min_consecutive_constraints-2-false") {
	{
		IP ip;
		vector<Sum> x;
		Sum x_sum = 0;
		for (int i = 1; i <= 3; ++i) {
			x.emplace_back(ip.add_boolean());
			x_sum += x.back();
		}
		ip.add_constraint(x_sum == 2);
		ip.add_min_consecutive_constraints(2, x, false);
		int num_solutions = 0;
		IPSolver solver;
		solver.set_silent(true);
		auto solutions = solver.solutions(&ip);
		while (solutions->get()) {
			num_solutions++;
		}
		REQUIRE(num_solutions == 2);
	}

	{
		IP ip;
		vector<Sum> x;
		Sum x_sum = 0;
		for (int i = 1; i <= 5; ++i) {
			x.emplace_back(ip.add_boolean());
			x_sum += x.back();
		}
		ip.add_constraint(x_sum == 4);
		ip.add_min_consecutive_constraints(2, x, false);
		int num_solutions = 0;
		IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
		auto solutions = solver.solutions(&ip);
		while (solutions->get()) {
			num_solutions++;
		}
		REQUIRE(num_solutions == 3);
	}
}

TEST_CASE("add_min_consecutive_constraints-2-true") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 3; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 2);
	ip.add_min_consecutive_constraints(2, x, true);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	REQUIRE(num_solutions == 3);
}

TEST_CASE("add_min_consecutive_constraints-3-false") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 5; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 3);
	ip.add_min_consecutive_constraints(3, x, false);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	CHECK(num_solutions == 3);
}

TEST_CASE("add_min_consecutive_constraints-3-true") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 5; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 3);
	ip.add_min_consecutive_constraints(3, x, true);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	REQUIRE(num_solutions == 5);
}

TEST_CASE("add_min_consecutive_constraints-5-false") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 6);
	ip.add_min_consecutive_constraints(5, x, false);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	CHECK(num_solutions == 7);
}

TEST_CASE("add_min_consecutive_constraints-5-true") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_constraint(x_sum == 3);
	ip.add_min_consecutive_constraints(5, x, true);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	REQUIRE(num_solutions == 4);
}

TEST_CASE("add_min_consecutive_constraints-all-true") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_min_consecutive_constraints(12, x, true);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		num_solutions++;
	}
	REQUIRE(num_solutions == 1);
}

TEST_CASE("add_min_max_consecutive_constraints-4") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_min_consecutive_constraints(4, x, false);
	ip.add_max_consecutive_constraints(4, x);

	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		CHECK(x_sum.value() <= 8);
		int consecutive = 0;
		for (auto& xx : x) {
			if (xx.value() > 0.5) {
				consecutive++;
			} else {
				consecutive = 0;
			}
			CHECK(consecutive <= 4);
		}
		num_solutions++;
	}
	CHECK(num_solutions == 20);
}

TEST_CASE("add_min_max_consecutive_constraints-5") {
	IP ip;
	vector<Sum> x;
	Sum x_sum = 0;
	for (int i = 1; i <= 12; ++i) {
		x.emplace_back(ip.add_boolean());
		x_sum += x.back();
	}
	ip.add_min_consecutive_constraints(5, x, false);
	ip.add_max_consecutive_constraints(5, x);
	int num_solutions = 0;
	IpToSatSolver solver(bind(minisat_solver, true));  // Faster.
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		int consecutive = 0;
		for (auto& xx : x) {
			if (xx.value() > 0.5) {
			} else {
				consecutive = 0;
			}
			CHECK(consecutive <= 5);
		}
		num_solutions++;
	}
	CHECK(num_solutions == 12);
}

TEST_CASE("solve_relaxation") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	Sum objective = -x - y;
	ip.add_objective(objective);
	ip.add_constraint(2 * x <= 1);
	ip.add_constraint(2 * y <= 1);

	REQUIRE(IPSolver().solve_relaxation(&ip));
	CHECK(Approx(objective.value()) == -1);
	CHECK(Approx(x.value()) == 0.5);
	CHECK(Approx(y.value()) == 0.5);
	REQUIRE(solve(&ip));
	CHECK(Approx(objective.value()) == 0);
	CHECK(Approx(x.value()) == 0);
	CHECK(Approx(y.value()) == 0);
}

TEST_CASE("trivial_constraint") {
	IP ip;
	Sum s = 0;
	CHECK_FALSE(ip.add_constraint(s <= 1).is_valid());
	CHECK_THROWS(ip.add_constraint(s <= -1));
	CHECK_THROWS(ip.add_constraint(-1, s, -1));
	CHECK_THROWS(ip.add_constraint(1, s, 1));
	auto x = ip.add_boolean();
	CHECK_FALSE(ip.add_constraint(x <= 1).is_valid());
	auto y = ip.add_boolean();
	CHECK(ip.add_constraint(x + y <= 1).is_valid());
}

TEST_CASE("constraint_always_false") {
	IP ip;
	CHECK_THROWS(ip.add_constraint(Sum(1) == 0));
	auto x = ip.add_boolean();
	CHECK_THROWS(ip.add_constraint(0 * x == 1));
	CHECK_FALSE(ip.add_constraint(x == 1).is_valid());
}

TEST_CASE("single_variable_constraint") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	ip.add_constraint(x >= 1);
	ip.add_objective(x);

	auto y = ip.add_variable(IP::Real);
	ip.add_constraint(-y <= 1);
	ip.add_objective(y);

	REQUIRE(solve(&ip));
	CHECK(x.value() == 1);
	CHECK(y.value() == -1);
}

TEST_CASE("sum_without_variables") {
	Sum s = 4;
	CHECK(s.value() == 4);
}

TEST_CASE("exists_basics") {
	IP ip;

	int count = 0;
	for (auto i : ip.exists({1, 2, 3})) {
		count++;
		CHECK(i == count);
	}
	CHECK(count == 3);

	std::vector<double> v = {4.5, 4.5, 4.5, 4.5};
	count = 0;
	for (auto x : ip.exists(v)) {
		CHECK(x == 4.5);
		count++;
	}
	CHECK(count == 4);
}

TEST_CASE("exists_reformulation_with_extra_var") {
	// Like the previous test case, but with an extra variable.
	IP ip;
	auto Digits = {0, 1, 2, 3, 4};
	auto d = ip.add_vector(5, IP::Integer);
	for (auto& digit : d) {
		ip.add_bounds(0, digit, 9);
	}
	auto original = ip.add_variable(IP::Real);

	Sum original_sum = 0;
	for (int i : Digits) {
		original_sum += d[i] * pow(10, i);
	}
	ip.add_constraint(original == original_sum);

	for (int i : ip.exists(Digits)) {
		Sum new_sum = 0;
		for (int j : Digits) {
			if (j != i) {
				int exp = j < i ? j : j - 1;
				new_sum += d[j] * pow(10, exp);
			}
		}
		ip.add_constraint(original + new_sum == 41751);
	}

	REQUIRE(solve(&ip));
	CHECK(original.value() == 37956);
}

TEST_CASE("exists_empty") {
	IP ip;
	auto x = ip.add_variable();
	ip.add_constraint(x == 0);
	for (auto y : ip.exists({"String1", "String2", "String3"})) {
		// Nothing.
	}
	REQUIRE(solve(&ip));
}

TEST_CASE("exists_infeasible") {
	IP ip;
	auto x = ip.add_variable();
	ip.add_constraint(x == 0);
	for (auto y : ip.exists({"String1", "String2", "String3"})) {
		ip.add_constraint(x == 1);
	}
	REQUIRE(!solve(&ip));
}

TEST_CASE("exists_identical") {
	IP ip;
	auto x = ip.add_variable();
	ip.add_constraint(x == 1);
	for (auto y : ip.exists({1, 2})) {
		ip.add_constraint(x == 1);
	}
	REQUIRE(solve(&ip));
}

TEST_CASE("exists_identical_no_other_constraints_0") {
	IP ip;
	auto x = ip.add_variable(IP::Boolean);
	for (auto y : ip.exists({1, 2})) {
		ip.add_constraint(x == 0);
	}
	REQUIRE(solve(&ip));
}

TEST_CASE("exists_identical_no_other_constraints_1") {
	IP ip;
	auto x = ip.add_variable(IP::Boolean);
	for (auto y : ip.exists({1, 2})) {
		ip.add_constraint(x == 1);
	}
	REQUIRE(solve(&ip));
}

TEST_CASE("exists_two_clauses") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	for (auto i : ip.exists({1, 2})) {
		ip.add_constraint(x == i);
	}
	for (auto i : ip.exists({2, 3})) {
		ip.add_constraint(x == i);
	}
	REQUIRE(solve(&ip));
	CHECK(Approx(x.value()) == 2);
}

TEST_CASE("exists_three_clauses_infeasible") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	for (auto i : ip.exists({1, 2})) {
		ip.add_constraint(x == i);
	}
	for (auto i : ip.exists({2, 3})) {
		ip.add_constraint(x == i);
	}
	for (auto i : ip.exists({3, 4})) {
		ip.add_constraint(x == i);
	}
	REQUIRE(!solve(&ip));
}

TEST_CASE("exists_multiple_constraints") {
	IP ip;
	auto x = ip.add_variable(IP::Boolean);
	auto y = ip.add_variable(IP::Boolean);
	for (auto i : ip.exists({0, 1})) {
		ip.add_constraint(x == i);
		ip.add_constraint(y == 1 - i);
	}
	REQUIRE(solve(&ip));
	CHECK(Approx(x.value() + y.value()) == 1);
}

TEST_CASE("exists_multiple_constraints_2") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	auto y = ip.add_variable(IP::Real);
	auto z = ip.add_variable(IP::Real);
	for (auto i : ip.exists({7})) {
		std::clog << "Petter:" << i << "\n";
		ip.add_constraint(x == i);
		ip.add_constraint(y == i);
		ip.add_constraint(z == i);
	}
	REQUIRE(solve(&ip));
	CHECK(Approx(x.value()) == 7);
	CHECK(Approx(y.value()) == 7);
	CHECK(Approx(z.value()) == 7);
}

TEST_CASE("exists_nested") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	for (auto i : ip.exists({7})) {
		ip.add_constraint(x == i);
		CHECK_THROWS([&]() {
			for (auto j : ip.exists({7})) {
				ip.add_constraint(x == j);
			}
		}());
	}
}

TEST_CASE("callback") {
	IP ip;
	create_soduku_IP(ip, 2);
	int counter = 0;
	IPSolver solver;
	solver.set_callback([&counter, &ip]() {
		minimum_core_assert(ip.is_feasible_and_integral());
		counter++;
	});
	REQUIRE(solver.solutions(&ip)->get());
	CHECK(counter >= 1);
}

TEST_CASE("max_objective") {
	{
		IP ip;
		auto x1 = ip.add_variable(IP::Real);
		auto x2 = ip.add_variable(IP::Real);
		auto y = ip.max(x1, x2);
		ip.add_objective(y);

		ip.add_bounds(4, x1, 8);
		ip.add_bounds(5, x2, 7);
		REQUIRE(solve(&ip));
		CHECK(y.value() == 5);
	}

	{
		IP ip;
		auto x1 = ip.add_variable(IP::Real);
		auto x2 = ip.add_variable(IP::Real);
		auto y = ip.max(x1, x2);
		ip.add_objective(y);

		ip.add_bounds(-8, x1, -4);
		ip.add_bounds(-10, x2, -4);
		REQUIRE(solve(&ip));
		CHECK(y.value() == -8);
	}
}

TEST_CASE("abs_objective") {
	for (auto bounds : {make_pair(4, 8), make_pair(-8, -4)}) {
		IP ip;
		auto x = ip.add_variable(IP::Real);
		auto y = ip.abs(x);
		ip.add_objective(y);

		ip.add_bounds(bounds.first, x, bounds.second);
		REQUIRE(solve(&ip));
		CHECK(y.value() == 4);
	}
}

TEST_CASE("convex_objective_fail") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	auto y = ip.abs(x);
	CHECK_NOTHROW(ip.add_objective(y));
	CHECK_THROWS(ip.add_objective(-y));
}

TEST_CASE("convex_constraint1") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	ip.add_constraint(ip.abs(x) <= 5);
	ip.add_objective(x);
	REQUIRE(solve(&ip));
	CHECK(x.value() == -5);
}

TEST_CASE("convex_constraint2") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	ip.add_constraint(5 >= ip.abs(x));
	ip.add_objective(x);
	REQUIRE(solve(&ip));
	CHECK(x.value() == -5);
}

TEST_CASE("convex_constraint_fail") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	CHECK_THROWS(ip.add_constraint(ip.abs(x) >= 5));
}

TEST_CASE("dual_variables") {
	IP ip;
	auto x = ip.add_variable(IP::Real);
	auto y = ip.add_variable(IP::Real);
	ip.add_objective(-x - 2 * y);
	ip.add_constraint(x >= 0);
	auto invalid_dual = ip.add_constraint(y >= 0);
	auto dual = ip.add_constraint(x + y <= 1);
	CHECK(dual.is_valid());
	CHECK(!invalid_dual.is_valid());
	CHECK_FALSE(dual.is_available());

	REQUIRE(solve(&ip));
	REQUIRE(dual.is_available());
	CAPTURE(dual.value());
	CHECK(Approx(abs(dual.value())) == 2);
}

TEST_CASE("warm-start") {
	IP ip1;
	auto x1 = ip1.add_variable(IP::Real);
	auto x2 = ip1.add_variable(IP::Real);
	auto y = ip1.add_variable(IP::Real);
	ip1.add_objective(y);

	ip1.add_bounds(4, x1, 8);
	ip1.add_bounds(5, x2, 7);
	ip1.add_constraint(x1 + x2 + y >= 4);

	IPSolver solver;
	auto solutions = solver.solutions(&ip1);
	REQUIRE(solutions->get());
	CHECK(y.value() == -11);

	{
		IP ip2;
		auto x1 = ip2.add_variable(IP::Real);
		auto x2 = ip2.add_variable(IP::Real);
		auto y = ip2.add_variable(IP::Real);
		ip2.add_objective(y);

		ip2.add_bounds(4, x1, 8);
		ip2.add_bounds(5, x2, 7);
		ip2.add_constraint(x1 + x2 + y >= 3);
		solutions->warm_start(&ip2);
		REQUIRE(solutions->get());
		CHECK(y.value() == -12);
	}
}
