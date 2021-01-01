// Petter Strandmark.

#include <catch.hpp>

#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/minisatp.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
#include <minimum/linear/testing.h>

using namespace std;
using namespace minimum::linear;

void REQUIRE_NUM_SOLUTIONS(IP& ip, const Solver& solver, int num) {
	auto solutions = solver.solutions(&ip);

	int num_solutions = 0;
	while (solutions->get()) {
		num_solutions++;
	}
	CHECK(num_solutions == num);
}

void REQUIRE_NUM_SOLUTIONS(IP& ip, function<unique_ptr<SatSolver>()> sat_solver_factory, int num) {
	IpToSatSolver solver(std::move(sat_solver_factory));
	REQUIRE_NUM_SOLUTIONS(ip, solver, num);
}

bool solve(IP* ip, function<unique_ptr<SatSolver>()> sat_solver_factory) {
	IpToSatSolver solver(sat_solver_factory);
	return solver.solutions(ip)->get();
}

void simple_test(function<unique_ptr<SatSolver>()> sat_solver_factory) {
	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y == 1);
		REQUIRE(solve(&ip, sat_solver_factory));
		CHECK((x.value() + y.value() == 1));

		auto z = ip.add_boolean();
		ip.add_constraint(x + z == 1);
		ip.add_constraint(y + z == 1);
		CHECK(!solve(&ip, sat_solver_factory));
		ip.add_constraint(x + y + z == 1);
		CHECK(!solve(&ip, sat_solver_factory));
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		ip.add_objective(-2 * x + 10);
		ip.add_constraint(2 * x >= 0);
		CHECK(solve(&ip, sat_solver_factory));
		CHECK(x.value() == 1);
	}

	{
		IP ip;
		IpToSatSolver solver(bind(minisat_solver, true));
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_objective(2 * x + y);
		ip.add_constraint(x + y <= 1);
		solver.allow_ignoring_cost_function = true;
		CHECK(solver.solutions(&ip)->get());
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y <= 1);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 3);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		auto u = ip.add_boolean();
		auto v = ip.add_boolean();
		ip.add_constraint(x + y + z + w + u + v <= 5);
		ip.add_constraint(x + y + z + w + u + v >= 2);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 56);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + y == 2);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 1);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x + 2 * y == 1);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 1);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_constraint(x - y == 1);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 1);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		ip.add_constraint(x + y - z - w == 0);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 6);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		ip.add_constraint(x + y + z == 1);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 3);
	}

	{
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		ip.add_bounds(0, x, 0);
		ip.add_bounds(1, y, 1);
		REQUIRE(solve(&ip, sat_solver_factory));
		CHECK(x.bool_value() == 0);
		CHECK(y.bool_value() == 1);
		REQUIRE_NUM_SOLUTIONS(ip, sat_solver_factory, 1);
	}

	{
		IP ip;
		IpToSatSolver solver(bind(minisat_solver, true));
		solver.set_silent(true);
		auto x = ip.add_boolean();
		ip.add_objective(0.5 * x);
		CHECK_THROWS(solver.solutions(&ip)->get());
		solver.allow_ignoring_cost_function = true;
		CHECK(solver.solutions(&ip)->get());
	}
}

TEST_CASE("minisat") {
	simple_test([]() { return minisat_solver(true); });
}

TEST_CASE("minisatp") {
	// Disabled because Minisatp leaks memory.
	if (false) {
		IP ip;
		auto x = ip.add_boolean();
		auto y = ip.add_boolean();
		auto z = ip.add_boolean();
		auto w = ip.add_boolean();
		auto u = ip.add_boolean();
		auto v = ip.add_boolean();
		ip.add_constraint(x + y + z + w + u + v <= 5);
		ip.add_constraint(x + y + z + w + u + v >= 2);

		REQUIRE(MinisatpSolver{}.solutions(&ip)->get());
		CHECK(ip.is_feasible_and_integral());
	}
}

TEST_CASE("glucose") { simple_test(glucose_solver); }

TEST_CASE("simple_integer_programming_tests_minisat") {
	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	// "square" takes more than a second.
	CHECK_NOTHROW(simple_integer_programming_tests(solver, {"unbounded", "exists", "square"}));
}

TEST_CASE("simple_integer_programming_tests_glucose") {
	IpToSatSolver solver(glucose_solver);
	solver.set_silent(true);
	// "square" takes more than a second.
	CHECK_NOTHROW(simple_integer_programming_tests(solver, {"unbounded", "exists", "square"}));
}

TEST_CASE("sudoku") {
	IP ip;
	int n = 3;
	auto P = create_soduku_IP(ip, n);

	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	REQUIRE(solver.solutions(&ip)->get());

	vector<vector<int>> solution(n * n);

	cout << endl;
	for (int i = 0; i < n * n; ++i) {
		for (int j = 0; j < n * n; ++j) {
			solution[i].emplace_back();

			for (int k = 0; k < n * n; ++k) {
				if (P[i][j][k].bool_value()) {
					solution[i][j] = k + 1;
				}
			}
		}
	}

	for (int i = 0; i < n * n; ++i) {
		int row_sum = 0;
		int col_sum = 0;
		for (int j = 0; j < n * n; ++j) {
			row_sum += solution[i][j];
			col_sum += solution[j][i];
		}
		CHECK(row_sum == 45);
		CHECK(col_sum == 45);
	}
}

TEST_CASE("infeasible_on_creation_1") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x + y >= 10);
	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	CHECK_FALSE(solver.solutions(&ip)->get());
}

TEST_CASE("infeasible_on_creation_2") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x + y <= -10);
	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	CHECK_FALSE(solver.solutions(&ip)->get());
}

TEST_CASE("repeated") {
	IP ip;
	auto x = ip.add_boolean();
	ip.add_boolean();
	ip.add_objective(-3 * x);
	REQUIRE(solve_minisat(&ip));
	CHECK(x.value() == 1);
	ip.add_objective(4 * x);
	REQUIRE(solve_minisat(&ip));
	CHECK(x.value() == 0);
}

TEST_CASE("minisat-large") {
	IP ip;
	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	solver.allow_ignoring_cost_function = true;
	Sum sum = 0;
	for (int i = 0; i < 100; ++i) {
		auto x = ip.add_boolean(1);
		sum += x;
	}
	ip.add_constraint(sum == 50);
	CHECK(solver.solutions(&ip)->get());
	CHECK(sum.value() == 50);
}

TEST_CASE("sat-objective") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	ip.add_objective(4 * x + y + 2 * z + w);
	ip.add_constraint(x + y + z + w == 2);
	CHECK(solve_minisat(&ip));
	CHECK(!x.bool_value());
	CHECK(y.bool_value());
	CHECK(!z.bool_value());
	CHECK(w.bool_value());
}

TEST_CASE("sat-constraints_coefs_larger_than_1") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(4 * x - y - z >= 0);
	REQUIRE(solve_minisat(&ip));

	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 5);
}

TEST_CASE("sat-constraints_coefs_zero") {
	IP ip;
	auto x = ip.add_boolean();
	ip.add_constraint(0 * x >= 0);

	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 2);
}

TEST_CASE("sat-constraints_fractional_coeff-1") {
	IP ip;
	auto x = ip.add_boolean();
	ip.add_constraint(0.5 * x == 0.25);

	CHECK_THROWS(solve_minisat(&ip));
}

TEST_CASE("sat-constraints_fractional_coeff-2") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(0.5 * x + y >= 0);

	CHECK_THROWS(solve_minisat(&ip));
}

TEST_CASE("sat-constraints_coefs_less_than_minus_1") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(8 * x - 2 * y - 2 * z >= 0);

	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 5);
}

TEST_CASE("sat-constraints_coefs_negative_multiple_times") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	ip.add_constraint(8 * x - y - y - z - z >= 0);

	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 5);
}

TEST_CASE("sat-negative-objective") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	ip.add_objective(-4 * x - y - 2 * z - w);
	ip.add_constraint(x + y + z + w == 2);
	CHECK(solve_minisat(&ip));
	CHECK(x.bool_value());
	CHECK(!y.bool_value());
	CHECK(z.bool_value());
	CHECK(!w.bool_value());
}

TEST_CASE("sat-objective-next_solution") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();
	auto u = ip.add_boolean();
	auto v = ip.add_boolean();
	auto obj = 4 * x + y + 2 * z + w + u + v;
	ip.add_objective(obj);
	ip.add_constraint(x + y + z + w + u + v == 2);

	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 6);
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

	// 1 1 0
	// 0 1 1
	// 1 0 1
	REQUIRE_NUM_SOLUTIONS(
	    ip, []() { return minisat_solver(true); }, 3);
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
	IpToSatSolver solver(bind(minisat_solver, true));
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
	IpToSatSolver solver(bind(minisat_solver, true));
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

TEST_CASE("test_objective_positive") {
	IP ip;
	auto x = ip.add_boolean_vector(5);
	Sum objective = x[0] + x[1] + x[2] + x[3] + x[4];
	ip.add_objective(objective);
	ip.add_constraint(x[0] + x[1] + x[2] >= 2);
	ip.add_constraint(x[2] + x[3] + x[4] >= 2);
	CHECK(solve_minisat(&ip));
	CHECK(objective.value() == 3);
}

TEST_CASE("test_objective_negative") {
	IP ip;
	auto x = ip.add_boolean_vector(5);
	Sum objective = -x[0] - x[1] - x[2] - x[3] - x[4];
	ip.add_objective(objective);
	ip.add_constraint(x[0] + x[1] + x[2] <= 1);
	ip.add_constraint(x[2] + x[3] + x[4] <= 1);
	CHECK(solve_minisat(&ip));
	CHECK(objective.value() == -2);
}

TEST_CASE("test_objective_unsat") {
	IP ip;
	auto x = ip.add_boolean_vector(5);
	Sum objective = -x[0] - x[1] - x[2] - x[3] - x[4];
	ip.add_objective(objective);
	ip.add_constraint(x[0] + x[1] + x[2] <= 1);
	ip.add_constraint(x[0] + x[1] + x[2] >= 2);
	ip.add_constraint(x[2] + x[3] + x[4] <= 1);
	CHECK(!solve_minisat(&ip));
}

TEST_CASE("read_cnf") {
	auto data = R"(
c generated by sgen3tocnf
c nLits = 612
c actual nVars = 97
c Command: sgen3tocnf 024 0 2 24 97 4 2 0 /dev/null
p cnf 97 204
-1 -9 -53 0
-1 -9 -85 0
-1 -53 -85 0
-9 -53 -85 0
-37 -41 -54 0
-37 -41 -89 0
-37 -54 -89 0
-41 -54 -89 0
-10 -21 -49 0
-10 -21 -57 0
-10 -49 -57 0
-21 -49 -57 0
-2 -38 -65 0
-2 -38 -81 0
-2 -65 -81 0
-38 -65 -81 0
-22 -61 -69 0
-22 -61 -82 0
-22 -69 -82 0
-61 -69 -82 0
-62 -73 -77 0
-62 -73 -93 0
-62 -77 -93 0
-73 -77 -93 0
-5 -58 -78 0
-5 -58 -83 0
-5 -78 -83 0
-58 -78 -83 0
-25 -33 -42 0
-25 -33 -50 0
-25 -42 -50 0
-33 -42 -50 0
-11 -34 -70 0
-11 -34 -94 0
-11 -70 -94 0
-34 -70 -94 0
-13 -35 -55 0
-13 -35 -63 0
-13 -55 -63 0
-35 -55 -63 0
-17 -29 -56 0
-17 -29 -84 0
-17 -56 -84 0
-29 -56 -84 0
-3 -14 -18 0
-3 -14 -26 0
-3 -18 -26 0
-14 -18 -26 0
-6 -15 -66 0
-6 -15 -74 0
-6 -66 -74 0
-15 -66 -74 0
-39 -45 -71 0
-39 -45 -75 0
-39 -71 -75 0
-45 -71 -75 0
-4 -43 -46 0
-4 -43 -59 0
-4 -46 -59 0
-43 -46 -59 0
-19 -60 -76 0
-19 -60 -90 0
-19 -76 -90 0
-60 -76 -90 0
-16 -47 -51 0
-16 -47 -79 0
-16 -51 -79 0
-47 -51 -79 0
-12 -20 -64 0
-12 -20 -67 0
-12 -64 -67 0
-20 -64 -67 0
-7 -23 -36 0
-7 -23 -91 0
-7 -36 -91 0
-23 -36 -91 0
-27 -48 -86 0
-27 -48 -92 0
-27 -86 -92 0
-48 -86 -92 0
-24 -30 -80 0
-24 -30 -87 0
-24 -80 -87 0
-30 -80 -87 0
-28 -31 -68 0
-28 -31 -95 0
-28 -68 -95 0
-31 -68 -95 0
-8 -44 -72 0
-8 -44 -88 0
-8 -72 -88 0
-44 -72 -88 0
-32 -40 -52 0
-32 -40 -96 0
-32 -40 -97 0
-32 -52 -96 0
-32 -52 -97 0
-32 -96 -97 0
-40 -52 -96 0
-40 -52 -97 0
-40 -96 -97 0
-52 -96 -97 0
1 2 3 0
1 2 4 0
1 3 4 0
2 3 4 0
5 6 7 0
5 6 8 0
5 6 97 0
5 7 8 0
5 7 97 0
5 8 97 0
6 7 8 0
6 7 97 0
6 8 97 0
7 8 97 0
9 10 11 0
9 10 12 0
9 11 12 0
10 11 12 0
13 14 15 0
13 14 16 0
13 15 16 0
14 15 16 0
17 18 19 0
17 18 20 0
17 19 20 0
18 19 20 0
21 22 23 0
21 22 24 0
21 23 24 0
22 23 24 0
25 26 27 0
25 26 28 0
25 27 28 0
26 27 28 0
29 30 31 0
29 30 32 0
29 31 32 0
30 31 32 0
33 34 35 0
33 34 36 0
33 35 36 0
34 35 36 0
37 38 39 0
37 38 40 0
37 39 40 0
38 39 40 0
41 42 43 0
41 42 44 0
41 43 44 0
42 43 44 0
45 46 47 0
45 46 48 0
45 47 48 0
46 47 48 0
49 50 51 0
49 50 52 0
49 51 52 0
50 51 52 0
53 54 55 0
53 54 56 0
53 55 56 0
54 55 56 0
57 58 59 0
57 58 60 0
57 59 60 0
58 59 60 0
61 62 63 0
61 62 64 0
61 63 64 0
62 63 64 0
65 66 67 0
65 66 68 0
65 67 68 0
66 67 68 0
69 70 71 0
69 70 72 0
69 71 72 0
70 71 72 0
73 74 75 0
73 74 76 0
73 75 76 0
74 75 76 0
77 78 79 0
77 78 80 0
77 79 80 0
78 79 80 0
81 82 83 0
81 82 84 0
81 83 84 0
82 83 84 0
85 86 87 0
85 86 88 0
85 87 88 0
86 87 88 0
89 90 91 0
89 90 92 0
89 91 92 0
90 91 92 0
93 94 95 0
93 94 96 0
93 95 96 0
94 95 96 0
)";

	stringstream sin(data);
	auto ip = read_CNF(sin);
	CHECK(ip->get_number_of_variables() == 97);
	CHECK(IPSolver().solve_relaxation(ip.get()));
}
