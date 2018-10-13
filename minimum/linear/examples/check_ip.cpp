//
// Runs an IP saved in a proto file.
//
// For crash debugging of Coin solver.
//

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/proto.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program(int num_args, char* args[]) {
	using namespace std;
	check(num_args == 3, "Need a IP file name and solution file name.");

	Timer t("Reading file");
	ifstream in(args[1], ios::binary);
	check(bool(in), "Failed opening file.");
	in.seekg(0, in.end);
	int length = in.tellg();
	in.seekg(0, in.beg);
	vector<char> buffer(length);
	in.read(buffer.data(), length);
	check(bool(in), "Failed reading file.");

	IP ip(buffer.data(), buffer.size());

	t.next("Reading solution");
	proto::Solution solution;
	ifstream solution_in(args[2], ios::binary);
	check(solution.ParseFromIstream(&solution_in), "Could not parse solution file");
	t.OK();

	check(solution.primal_size() == ip.get_number_of_variables(),
	      "Solution primal size does not match.");
	check(solution.dual_size() == ip.get().constraint_size(), "Solution dual size does not match.");

	for (int i = 0; i < ip.get_number_of_variables(); ++i) {
		ip.set_solution(i, solution.primal(i));
	}
	for (int j = 0; j < ip.get().constraint_size(); ++j) {
		ip.set_dual_solution(j, solution.dual(j));
	}

	cout << "Number of primal variables: " << ip.get_number_of_variables() << "\n";
	cout << "Number of constraints: " << ip.get().constraint_size() << "\n";
	cout << "Current objective value: " << ip.get_entire_objective() << "\n";

	for (double eps = 1.0; eps >= 1e-16; eps *= 0.1) {
		bool feasible = ip.is_feasible(eps);
		if (!feasible) {
			cout << "NOT ";
		}
		cout << "feasible at 10^" << std::log10(eps) << "\n";
		if (!feasible) {
			break;
		}
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
