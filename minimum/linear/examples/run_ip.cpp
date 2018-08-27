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
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;
using namespace minimum::linear;

DEFINE_int32(n, 1, "Number of times to solve it.");

int main_program(int num_args, char* args[]) {
	using namespace std;

	Timer t("Reading file");
	check(num_args == 2, "Need a file name.");
	ifstream in(args[1], ios::binary);
	check(bool(in), "Failed opening file.");
	in.seekg(0, in.end);
	int length = in.tellg();
	in.seekg(0, in.beg);
	vector<char> buffer(length);
	in.read(buffer.data(), length);
	check(bool(in), "Failed reading file.");

	IP ip(buffer.data(), buffer.size());
	t.OK();

	t.next("Solving");
	IPSolver solver;
	solver.set_silent(true);
	for (int i = 0; i < FLAGS_n; ++i) {
		check(solver.solutions(&ip)->get(), "Solving failed.");
	}
	t.OK();

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
