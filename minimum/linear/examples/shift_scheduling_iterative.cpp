// Petter Strandmark 2014.
//
// Reads and solves problem files from “Employee Shift Scheduling Benchmark
// Data Sets”
// http://www.cs.nott.ac.uk/~tec/NRP/
//
// First argument to the program is the problem .txt file.
//
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/sqlite.h>
#include <minimum/core/string.h>
#include <minimum/core/time.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/scheduling_util.h>
using namespace minimum::core;
using namespace minimum::linear;

DEFINE_string(table, "solution", "Queries solution from this table. Default: solution.");
DEFINE_int32(start_day, 0, "Starts the refinement on this day.");
DEFINE_int32(start_solution, -1, "Starts the local solver from this solution.");

std::mt19937_64 rng;

int create_solve_and_update_solution(const proto::SchedulingProblem& problem,
                                     vector<vector<vector<int>>>& current_solution,
                                     const vector<vector<vector<Sum>>>& x,
                                     IP& ip) {
	int current_objective = objective_value(problem, current_solution);
	auto objective = create_ip(ip, problem, x);

	auto parse_solution = [&]() {
		for (int p = 0; p < problem.worker_size(); ++p) {
			for (int d = 0; d < problem.num_days(); ++d) {
				for (int s = 0; s < problem.shift_size(); ++s) {
					current_solution[p][d][s] = x[p][d][s].value() > 0.5 ? 1 : 0;
				}
			}
		}
	};

	IPSolver solver;
	solver.time_limit_in_seconds = 10.0;
	solver.set_silent(true);

	if (solver.solutions(&ip)->get()) {
		parse_solution();

		current_objective = int(objective.value() + 0.5);
		clog << "-- IP cost is " << current_objective << endl;
	} else {
		clog << "-- IP FAILED. Time-limit?" << endl;
	}

	return current_objective;
}

int optimize_single_window(const proto::SchedulingProblem& problem,
                           vector<vector<vector<int>>>& current_solution,
                           int start_day,
                           int end_day,
                           const vector<int>& current_staff_indices) {
	clog << "-- Optimizing for " << current_staff_indices.size() << " staff members." << endl;
	IP ip;
	vector<vector<vector<Sum>>> x;
	for (int p = 0; p < problem.worker_size(); ++p) {
		x.emplace_back();
		for (int d = 0; d < problem.num_days(); ++d) {
			x.back().emplace_back();

			for (int s = 0; s < problem.shift_size(); ++s) {
				if (find(begin(current_staff_indices), end(current_staff_indices), p)
				        != end(current_staff_indices)
				    && start_day <= d && d < end_day
				    && problem.worker(p).shift_limit(s).max() > 0) {
					x.back().back().emplace_back(ip.add_boolean());
				} else {
					x.back().back().emplace_back(current_solution[p][d][s]);
				}
			}
		}
	}

	return create_solve_and_update_solution(problem, current_solution, x, ip);
}

int perform_sliding_window_optimization(const proto::SchedulingProblem& problem,
                                        vector<vector<vector<int>>>& current_solution,
                                        int window_size = 6,
                                        int max_num_staff = 32,
                                        std::function<void()> print_function = nullptr) {
	window_size = min(window_size, problem.num_days());

	int objective_value = -1;

	static bool first_time = true;
	int first_day = 0;
	if (first_time) {
		first_day = FLAGS_start_day;
	}
	first_time = false;

	for (int start_day = first_day; start_day < problem.num_days() - window_size + 1;
	     start_day += 1) {
		clog << "Optimizing window [" << start_day << ", " << start_day + window_size << ") out of "
		     << problem.num_days() << " days." << endl;

		std::vector<int> staff_indices;
		for (int p = 0; p < problem.worker_size(); ++p) {
			staff_indices.push_back(p);
		}
		// shuffle(begin(staff_indices), end(staff_indices), rng);

		for (size_t start_index = 0; start_index < staff_indices.size();
		     start_index += max_num_staff) {
			size_t end_index = min(staff_indices.size(), start_index + max_num_staff);
			vector<int> current_staff_indices(begin(staff_indices) + start_index,
			                                  begin(staff_indices) + end_index);
			sort(begin(current_staff_indices), end(current_staff_indices));

			auto this_objective_value = optimize_single_window(problem,
			                                                   current_solution,
			                                                   start_day,
			                                                   start_day + window_size,
			                                                   current_staff_indices);
			if (this_objective_value < objective_value && print_function) {
				print_function();
			}
			objective_value = this_objective_value;
		}
	}
	return objective_value;
}

int main_program(int num_args, char* args[]) {
	using namespace std;

	minimum_core_assert(num_args >= 2);
	ifstream fin(args[1]);
	auto problem = read_nottingham_instance(fin);

	auto db = SqliteDb::fromFile(to_string(args[1], ".db"));
	db.make_statement(
	      "CREATE TABLE IF NOT EXISTS "
	      "refined ("
	      "objective INTEGER UNIQUE,"
	      "start_objective INTEGER,"
	      "time FLOAT,"
	      "solution TEXT,"
	      "options TEXT"
	      ");")
	    .execute();
	auto insert_solution =
	    db.make_statement("INSERT OR IGNORE INTO refined VALUES(?1, ?2, ?3, ?4, ?5);");

	double start_time = minimum::core::wall_time();

	Timer t("Reading initial solution from DB");
	std::string initial_solution;
	if (FLAGS_start_solution < 0) {
		initial_solution = db.make_statement<string>(to_string("SELECT solution FROM ",
		                                                       FLAGS_table,
		                                                       " ORDER BY objective ASC LIMIT 1;"))
		                       .execute()
		                       .get<0>();
	} else {
		initial_solution =
		    db.make_statement<string>(
		          to_string("SELECT solution FROM ", FLAGS_table, " WHERE objective = ?1 LIMIT 1;"))
		        .execute(FLAGS_start_solution)
		        .get<0>();
	}
	t.next("Loading solution");
	istringstream solution_file(initial_solution);
	auto current_solution = load_solution(solution_file, problem);
	t.next("Computing objective");
	auto start_objective = objective_value(problem, current_solution);
	t.OK();
	clog << "-- Starting solution is " << start_objective << "\n";

	auto print_solution = [&]() {
		double elapsed_time = minimum::core::wall_time() - start_time;
		auto obj = objective_value(problem, current_solution);

		clog << "-- Objective value is " << obj << endl;

		Timer t("Writing output to DB");
		stringstream sout;
		save_solution(sout, problem, current_solution);
		insert_solution.execute(obj, start_objective, elapsed_time, sout.str(), "");
		t.OK();
	};

	int window_size = 8;

	int max_num_staff = 16;
	if (problem.shift_size() >= 16) {
		max_num_staff = 5;
	}

	for (int iter = 1; iter <= 100000; ++iter) {
		Timer t("Trying to improve solution quickly");
		auto result = quick_solution_improvement(problem, current_solution);
		t.OK();
		if (result.empty()) {
			break;
		}
		clog << "-- " << result << endl;
	}

	print_solution();

	auto previous_best = numeric_limits<int>::max();
	for (int iter = 1; iter <= 1000; ++iter) {
		int curr = perform_sliding_window_optimization(
		    problem, current_solution, window_size, max_num_staff, print_solution);

		minimum_core_assert(curr <= previous_best);
		if (curr == previous_best) {
			break;
		}
		previous_best = curr;

		print_solution();
	}
	print_solution();

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
