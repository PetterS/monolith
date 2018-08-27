#include <algorithm>
using namespace std;

#include <catch.hpp>

#include <minimum/core/range.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/symmetry.h>
using namespace minimum::core;
using namespace minimum::linear;

template <template <typename, typename> class Container, typename T, typename Alloc>
bool contains(const Container<T, Alloc>& container, T arg) {
	return std::find(container.begin(), container.end(), arg) != container.end();
}

TEST_CASE("simple-6") {
	IP ip;
	auto x = ip.add_boolean_vector(4);
	ip.add_constraint(x[0] + x[1] == 1);
	ip.add_constraint(x[0] + x[2] == 1);
	ip.add_constraint(x[0] + x[3] == 1);
	auto orbits = analyze_symmetry(ip).orbits;
	CHECK(contains(orbits, {0}));
	CHECK(contains(orbits, {1, 2, 3}));
}

TEST_CASE("simple-2") {
	IP ip;
	auto x = ip.add_boolean_vector(4);
	ip.add_constraint(x[0] + x[1] == 1);
	ip.add_constraint(x[0] + x[2] == 1);
	ip.add_constraint(x[0] + x[3] == 1);
	ip.add_constraint(x[1] + x[2] == 5);
	auto orbits = analyze_symmetry(ip).orbits;
	CHECK(contains(orbits, {0}));
	CHECK(contains(orbits, {3}));
	CHECK(contains(orbits, {1, 2}));
}

TEST_CASE("square") {
	// 0, 1, 2
	// 3, 4, 5
	// 6, 7, 8
	//
	// 9
	const int n = 3;

	IP ip;
	auto value = ip.add_grid(n, n, IP::Integer);
	auto s = ip.add_variable(IP::Integer);
	ip.add_bounds(1, s, 100);
	for (int j : range(n)) {
		for (int i : range(n)) {
			ip.add_bounds(1, value[i][j], n * n);
		}
	}

	for (int i : range(n)) {
		Sum i_sum = 0;
		for (int j : range(n)) {
			i_sum += value[i][j];
		}
		ip.add_constraint(i_sum == s);
	}

	for (int j : range(n)) {
		Sum j_sum = 0;
		for (int i : range(n)) {
			j_sum += value[i][j];
		}
		ip.add_constraint(j_sum == s);
	}

	Sum diag1 = 0;
	Sum diag2 = 0;
	for (int i : range(n)) {
		diag1 += value[i][i];
		diag2 += value[n - i - 1][i];
	}
	ip.add_constraint(diag1 == s);
	ip.add_constraint(diag2 == s);

	auto orbits = analyze_symmetry(ip).orbits;
	CHECK(contains(orbits, {4}));
	CHECK(contains(orbits, {9}));
	CHECK(contains(orbits, {0, 2, 6, 8}));
	CHECK(contains(orbits, {1, 3, 5, 7}));
}

//     0   1   2
//   11  12  13  3
// 10  17  18  14  4
//   9  16  15   5
//     8   7   6
TEST_CASE("hexagon-integers") {
	IP ip;
	auto x = ip.add_vector(19, IP::Integer);
	for (int i : range(19)) {
		ip.add_bounds(1, x[i], 19);
	}
	auto rhs = 38;
	ip.add_constraint(x[0] + x[1] + x[2] == rhs);
	ip.add_constraint(x[11] + x[12] + x[13] + x[3] == rhs);
	ip.add_constraint(x[10] + x[17] + x[18] + x[14] + x[4] == rhs);
	ip.add_constraint(x[9] + x[16] + x[15] + x[5] == rhs);
	ip.add_constraint(x[8] + x[7] + x[6] == rhs);

	ip.add_constraint(x[10] + x[9] + x[8] == rhs);
	ip.add_constraint(x[11] + x[17] + x[16] + x[7] == rhs);
	ip.add_constraint(x[0] + x[12] + x[18] + x[15] + x[6] == rhs);
	ip.add_constraint(x[1] + x[13] + x[14] + x[5] == rhs);
	ip.add_constraint(x[2] + x[3] + x[4] == rhs);

	ip.add_constraint(x[0] + x[11] + x[10] == rhs);
	ip.add_constraint(x[1] + x[12] + x[17] + x[9] == rhs);
	ip.add_constraint(x[2] + x[13] + x[18] + x[16] + x[8] == rhs);
	ip.add_constraint(x[3] + x[14] + x[15] + x[7] == rhs);
	ip.add_constraint(x[4] + x[5] + x[6] == rhs);

	auto orbits = analyze_symmetry(ip).orbits;
	CHECK(contains(orbits, {18}));
	CHECK(contains(orbits, {0, 2, 4, 6, 8, 10}));
}

void no_test_case() {
	const int n = 3;

	IP ip;
	auto x = ip.add_boolean_cube(n, n, n * n);
	// auto s = ip.add_variable(IP::Integer);
	// ip.add_bounds(1, s, 1000);
	int s = ((n * n) * (n * n + 1)) / (2 * n);

	// Every square contains exactly one number.
	for (int i : range(n)) {
		for (int j : range(n)) {
			Sum k_sum = 0;
			for (int k : range(n * n)) {
				k_sum += x[i][j][k];
			}
			ip.add_constraint(k_sum == 1);
		}
	}

	// Every number occurs exactly once.
	for (int k : range(n * n)) {
		Sum ij_sum = 0;
		for (int i : range(n)) {
			for (int j : range(n)) {
				ij_sum += x[i][j][k];
			}
		}
		ip.add_constraint(ij_sum == 1);
	}

	// Create a helper variable containing the integer at (i, j).
	vector<vector<Sum>> value(n, vector<Sum>(n, 0));
	for (int i : range(n)) {
		for (int j : range(n)) {
			for (int k : range(n * n)) {
				value[i][j] += x[i][j][k] * (k + 1);
			}
		}
	}

	for (int i : range(n)) {
		Sum i_sum = 0;
		for (int j : range(n)) {
			i_sum += value[i][j];
		}
		ip.add_constraint(i_sum == s);
	}

	for (int j : range(n)) {
		Sum j_sum = 0;
		for (int i : range(n)) {
			j_sum += value[i][j];
		}
		ip.add_constraint(j_sum == s);
	}

	Sum diag1 = 0;
	Sum diag2 = 0;
	for (int i : range(n)) {
		diag1 += value[i][i];
		diag2 += value[n - i - 1][i];
	}
	ip.add_constraint(diag1 == s);
	ip.add_constraint(diag2 == s);

	auto orbits = analyze_symmetry(ip).orbits;
	if (n == 3) {
		for (int k : range(n * n)) {
			CHECK(contains(orbits, {4 * n * n + k}));
		}
	}
	remove_symmetries(&ip);

	IpToSatSolver solver(bind(minisat_solver, true));
	auto solutions = solver.solutions(&ip);

	int num_solutions = 0;
	while (solutions->get()) {
		num_solutions++;
	}
	CHECK(num_solutions == 1);
}
