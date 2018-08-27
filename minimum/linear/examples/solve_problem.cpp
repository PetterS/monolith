
#include <fstream>
#include <iostream>

#include <gflags/gflags.h>

#include <minimum/core/color.h>
#include <minimum/linear/solver_from_command_line.h>
using namespace minimum::linear;

int main_program(int argc, char* argv[]) {
	using namespace std;
	if (argc != 2) {
		cerr << "Need one file name.\n";
		return 2;
	}

	minimum::core::Timer t("Reading input file");
	ifstream fin(argv[1]);
	auto ip = read_MPS(fin);
	t.OK();

	t.next("Creating solver");
	auto solver = create_solver_from_command_line();
	t.OK();

	if (!solver->solutions(ip.get())->get()) {
		cerr << "Solver failed.\n";
		return 3;
	}

	return 0;
}

int main(int argc, char* argv[]) {
	try {
		std::string usage = "Solves a problem in MPS format.\n";
		usage += argv[0];
		usage += " <flags> <file name>";
		gflags::SetUsageMessage(usage);
		gflags::ParseCommandLineFlags(&argc, &argv, true);
		return main_program(argc, argv);
	} catch (std::exception& err) {
		std::cerr << "ERROR: " << err.what() << std::endl;
		return 1;
	}
}
