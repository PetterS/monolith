// Petter Strandmark.
//
// Reads and solves problem files from NSPLib
// http://www.projectmanagement.ugent.be/?q=research/personnel_scheduling/nsp
//
// First argument to the program is the problem .nsp file. The second
// argument is the .gen-file with the constraints to use.
//
#include <fstream>
#include <string>
#include <vector>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/grid.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/scheduling_util.h>
using namespace minimum::core;
using namespace minimum::linear;

int main_program(int num_args, char* args[]) {
	using namespace std;

	minimum_core_assert(num_args >= 3);
	fstream problem_data(args[1]);
	fstream case_data(args[2]);
	auto problem = read_gent_instance(problem_data, case_data);
	string solution_file_name = "";
	if (num_args >= 4) {
		solution_file_name = args[3];
	}

	string problem_number = args[1];
	auto p1 = problem_number.find_last_of("/\\");
	if (p1 != string::npos) {
		problem_number = problem_number.substr(p1 + 1);
	}
	auto p2 = problem_number.find(".");
	if (p2 != string::npos) {
		problem_number = problem_number.substr(0, p2);
	}

	clog << "Number of nurses: " << problem.worker_size() << endl;
	clog << "Number of days  : " << problem.num_days() << endl;
	clog << "Number of shifts: " << problem.shift_size() << endl;

	double start_time = minimum::core::wall_time();

	IP ip;

	// The assignment variables.
	//   x[n][d][s] = 1 ⇔ Nurse n is working day d, shift s.
	auto x = ip.add_boolean_cube(problem.worker_size(), problem.num_days(), problem.shift_size());
	create_ip(ip, problem, x);
	double startup_time = minimum::core::wall_time() - start_time;

	clog << "IP created." << endl;

	auto solution = make_grid<int>(problem.worker_size(), problem.num_days(), problem.shift_size());
	auto assign_solution = [&x, &problem, &solution]() {
		for (int n = 0; n < problem.worker_size(); ++n) {
			for (int d = 0; d < problem.num_days(); ++d) {
				for (int s = 0; s < problem.shift_size(); ++s) {
					solution[n][d][s] = x[n][d][s].bool_value() ? 1 : 0;
				}
			}
		}
	};

	start_time = minimum::core::wall_time();

	// Try SAT solver first.
	bool try_SAT_solver = false;
	if (try_SAT_solver) {
		start_time = minimum::core::wall_time();
		IpToSatSolver solver(bind(minisat_solver, false));
		solver.allow_ignoring_cost_function = true;

		auto solutions = solver.solutions(&ip);
		int num_solutions = 0;
		while (solutions->get()) {
			double current_time = minimum::core::wall_time() - start_time;
			num_solutions++;
			assign_solution();
			clog << "SAT solver found solution of " << objective_value(problem, solution) << " in "
			     << current_time << " seconds." << endl;
			if (num_solutions >= 1) {  // Increase for more solutions.
				break;
			}
		}
	}

	// Solve in completely silent mode to keep stdout clean.
	IPSolver solver;
	solver.time_limit_in_seconds = 5.0;
	solver.set_silent(true);

	if (solver.solutions(&ip)->get()) {
		clog << "IP solved." << endl;
		assign_solution();
	} else if (!try_SAT_solver) {
		clog << "IP was not solved." << endl;
		cout << problem_number << "\t"
		     << "infeasible/timelimit"
		     << "\t" << 0 << endl;
		return 2;
	}

	double solver_time = minimum::core::wall_time() - start_time;

	if (solution_file_name != "") {
		ofstream fout(solution_file_name);
		save_solution(fout, problem, solution);
	} else {
		print_solution(cout, problem, solution);
	}

	clog << "Obtained objective function " << objective_value(problem, solution) << " in "
	     << solver_time << " seconds." << endl;

	cout << problem_number << "\t" << objective_value(problem, solution) << "\t"
	     << startup_time + solver_time << endl;

	clog << endl;
	return 0;
}

int main(int num_args, char* args[]) {
	gflags::ParseCommandLineFlags(&num_args, &args, true);
	try {
		return main_program(num_args, args);
	} catch (std::exception& err) {
		std::cerr << "Error: " << err.what() << std::endl;
		return 1;
	}
}
