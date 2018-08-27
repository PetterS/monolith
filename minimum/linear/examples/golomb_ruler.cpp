//
// A solver for Golomb rulers.
// http://en.wikipedia.org/wiki/Golomb_ruler
//

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/constraint_solver.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;
	IP ip;
	double start_time = wall_time();

	int max_length = 35;         // Maximum length of the ruler.
	int order = 8;               // The order of the created ruler.
	bool optimize_order = true;  // If true, instead maximize the order for the specific length.

	// Whether there is a mark at a certain position.
	int num_vars = max_length + 1;
	auto mark = ip.add_boolean_vector(num_vars);
	// Require a mark at position 0.
	ip.add_bounds(1, mark[0], 1);

	// Either maximize the order or find a ruler with a specific one.
	if (optimize_order) {
		ip.add_objective(-sum(mark));
	} else {
		ip.add_constraint(sum(mark) == order);
	}

	// For each distance, require that it occurs at most once.
	for (int dist = 1; dist < num_vars; ++dist) {
		PseudoBoolean distance_occurrences;
		for (int i : range(num_vars)) {
			for (int j : range(i + 1, num_vars)) {
				if (j - i == dist) {
					distance_occurrences += mark[i] * mark[j];
				}
			}
		}
		ip.add_pseudoboolean_constraint(distance_occurrences <= 1);
	}

	double elapsed_time = wall_time() - start_time;
	cout << "Building model: " << elapsed_time << "s.\n";
	start_time = wall_time();

	IpToSatSolver solver(bind(minisat_solver, true));
	solver.set_silent(true);
	// ConstraintSolver solver;

	auto solutions = solver.solutions(&ip);
	check(solutions->get(), "Could not solve problem.");
	elapsed_time = wall_time() - start_time;
	cout << "Solving model: " << elapsed_time << "s.\n";
	do {
		for (int i = 0; i < num_vars; ++i) {
			if (mark[i].bool_value()) {
				cout << i << " ";
			}
		}
		if (optimize_order) {
			cout << "order=" << sum(mark).value() << " ";
		}
		cout << "\n";
	} while (solutions->get());

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
