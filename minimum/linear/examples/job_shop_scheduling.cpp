// Description from GLPK:
//
// The Job-Shop Scheduling Problem (JSSP) is to schedule a set of jobs
// on a set of machines, subject to the constraint that each machine can
// handle at most one job at a time and the fact that each job has a
// specified processing order through the machines. The objective is to
// schedule the jobs so as to minimize the maximum of their completion
// times.
//
// Reference:
// D. Applegate and W. Cook, "A Computational Study of the Job-Shop
// Scheduling Problem", ORSA J. On Comput., Vol. 3, No. 2, Spring 1991,
// pp. 149-156.

#include <iomanip>
#include <iostream>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program() {
	using namespace std;

	// Permutation of the machines, which represents the processing order
	// of j through the machines: j must be processed first on sigma[j,1],
	// then on sigma[j,2], etc.
	vector<vector<int>> sigma = {{3, 1, 2, 4, 6, 5},
	                             {2, 3, 5, 6, 1, 4},
	                             {3, 4, 6, 1, 2, 5},
	                             {2, 1, 3, 4, 5, 6},
	                             {3, 2, 5, 6, 1, 4},
	                             {2, 4, 6, 1, 5, 3}};

	// Processing time of j on a.
	vector<vector<int>> p = {{3, 6, 1, 7, 6, 3},
	                         {10, 8, 5, 4, 10, 10},
	                         {9, 1, 5, 4, 7, 8},
	                         {5, 5, 5, 3, 8, 9},
	                         {3, 3, 9, 1, 5, 4},
	                         {10, 3, 1, 3, 4, 9}};

	int num_jobs = sigma.size();
	int num_machines = sigma[0].size();

	IP ip;

	// Starting time of j on a.
	auto x = ip.add_grid(num_jobs, num_machines, IP::Real);
	for (int j : range(num_jobs)) {
		for (int a : range(num_machines)) {
			ip.add_constraint(x[j][a] >= 0);
			if (a > 0) {
				// j can be processed on sigma[j,a] only after it has been completely
				// processed on sigma[j, a - 1]
				ip.add_constraint(x[j][sigma[j][a] - 1]
				                  >= x[j][sigma[j][a - 1] - 1] + p[j][sigma[j][a - 1] - 1]);
			}
		}
	}

	// The disjunctive condition that each machine can handle at most one
	// job at a time is the following:
	//
	// x[i,a] >= x[j,a] + p[j,a]  or  x[j,a] >= x[i,a] + p[i,a]
	//
	// for all i, j in J, a in M. This condition is modeled through binary
	// variables Y as shown below.
	auto Y = ip.add_boolean_cube(num_jobs, num_jobs, num_machines);
	int K = 0;
	for (auto& row : p) {
		for (auto val : row) {
			K += val;
		}
	}
	for (int j : range(num_jobs)) {
		for (int i : range(num_jobs)) {
			if (i != j) {
				for (int a : range(num_machines)) {
					ip.add_constraint(x[i][a] >= x[j][a] + p[j][a] - K * Y[i][j][a]);
					ip.add_constraint(x[j][a] >= x[i][a] + p[i][a] - K * (1 - Y[i][j][a]));
				}
			}
		}
	}

	auto z = ip.add_variable(IP::Real);
	for (int j : range(num_jobs)) {
		ip.add_constraint(z >= x[j][sigma[j][num_machines - 1] - 1]
		                           + p[j][sigma[j][num_machines - 1] - 1]);
	}
	ip.add_objective(z);
	check(solve(&ip), "Could not solve problem.");
	cout << z.value() << " is the completion time.\n";
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
