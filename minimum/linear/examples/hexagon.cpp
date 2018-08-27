// A model that is easy to solve with constraint programming
// but difficult with open source IP solvers.
//
// Compute a magic hexagon with these indices:
//
//     0   1   2
//   11  12  13  3
// 10  17  18  14  4
//   9  16  15   5
//     8   7   6
//
// With Cbc:
// Total time (CPU seconds):       166.50   (Wallclock seconds):       166.50
//
// Minisatp works best here.
#include <iostream>
#include <numeric>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/constraint_solver.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;
	IP ip;

	vector<Sum> x(19, 0);
	vector<Sum> number(19, 0);
	for (int i : range(19)) {
		auto xs = ip.add_boolean_vector(19);
		for (int j : range(19)) {
			x[i] += (j + 1) * xs[j];
			number[j] += xs[j];
		}
		ip.add_constraint(accumulate(xs.begin(), xs.end(), Sum(0)) == 1);
	}

	for (int i : range(19)) {
		ip.add_constraint(number[i] == 1);
	}

	auto rhs = 38;

	// auto rhs = ip.add_variable(IP::Integer);
	// ip.add_bounds(1, rhs, 100);
	// ip.add_constraint(x[2] == 3);
	// ip.add_constraint(x[3] == 17);

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

	ConstraintSolver solver;
	auto solutions = solver.solutions(&ip);
	while (solutions->get()) {
		cout << "\nSolution:\n";
		printf("    %2.0f  %2.0f  %2.0f\n", x[0].value(), x[1].value(), x[2].value());
		printf("  %2.0f  %2.0f  %2.0f  %2.0f\n",
		       x[11].value(),
		       x[12].value(),
		       x[13].value(),
		       x[3].value());
		printf("%2.0f  %2.0f  %2.0f  %2.0f  %2.0f\n",
		       x[10].value(),
		       x[17].value(),
		       x[18].value(),
		       x[14].value(),
		       x[4].value());
		printf("  %2.0f  %2.0f  %2.0f  %2.0f\n",
		       x[9].value(),
		       x[16].value(),
		       x[15].value(),
		       x[5].value());
		printf("    %2.0f  %2.0f  %2.0f\n", x[8].value(), x[7].value(), x[6].value());
	}

	return 0;
}

int main() {
	try {
		return main_program();
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
