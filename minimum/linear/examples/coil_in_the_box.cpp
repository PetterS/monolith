// Solves the coil-in-the-box problem.
// See https://en.wikipedia.org/wiki/Snake-in-the-box
#include <cstdlib>
#include <iostream>

#include <gflags/gflags.h>

#include <minimum/core/main.h>
#include <minimum/linear/connected_component.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/solver.h>
using namespace minimum::linear;

DEFINE_int32(n, 5, "Dimension of hypercube.");
DEFINE_string(mps_file_name, "", "If set, saves the IP as an MPS file.");

int main_program(int, char* []) {
	const int n = FLAGS_n;
	if (n < 2) {
		return 2;
	}
	IP ip;
	const int npoints = 1 << n;
	auto x = ip.add_boolean_vector(npoints, -1);

	for (int p = 0; p < npoints; ++p) {
		Sum neighbor_sum = 0;
		for (int i = 0; i < n; ++i) {
			unsigned neighbor = p ^ (1 << i);
			neighbor_sum += x.at(neighbor);
		}

		ip.add_constraint(neighbor_sum >= 2 - n * (1 - x.at(p)));
		ip.add_constraint(neighbor_sum <= 2 + n * (1 - x.at(p)));
	}
	ip.add_constraint(x.at(0) == 1);
	ip.add_constraint(x.at(1) == 1);
	ip.add_constraint(x.at(3) == 1);

	ConnectedComponentFlow connectivity(&ip);
	// ConnectedComponentNoFlow connectivity(&ip, 30);
	for (int p = 0; p < npoints; ++p) {
		connectivity.add_node(x.at(p));
	}
	for (int p = 0; p < npoints; ++p) {
		for (int i = 0; i < n; ++i) {
			unsigned neighbor = p ^ (1 << i);
			connectivity.add_edge(x.at(p), x.at(neighbor));
		}
	}
	connectivity.set_first_node(x.at(0));
	connectivity.finish();

	IPSolver solver;
	if (!FLAGS_mps_file_name.empty()) {
		solver.save_MPS(ip, FLAGS_mps_file_name);
	}

	// IpToSatSolver solver([]() { return minisat_solver(true); });
	// solver.set_silent(true);

	auto solutions = solver.solutions(&ip);
	if (solutions->get()) {
		int L = 0;
		for (int p = 0; p < npoints; ++p) {
			L += x[p].bool_value() ? 1 : 0;
		}
		std::cout << "n=" << n << " L=" << L << "\n";
	} else {
		std::cout << "Error. Did not find a solution.";
		return 1;
	}
	return 0;
}

int main(int num_args, char* args[]) { main_runner(main_program, num_args, args); }
