#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <string_view>

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/main.h>
#include <minimum/core/random.h>
#include <minimum/core/sqlite.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/retail_scheduling_pricing.h>
#include <minimum/linear/retail_scheduling.h>
#include <minimum/linear/scheduling_util.h>
#include <version/version.h>

using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;
using namespace std;

string get_all_command_line_options() {
	stringstream sout_flags;
	vector<gflags::CommandLineFlagInfo> all_flags;
	gflags::GetAllFlags(&all_flags);
	for (auto option : all_flags) {
		if (!option.is_default) {
			sout_flags << option.name << "=" << option.current_value << " ";
		}
	}
	return sout_flags.str();
}

int main_program(int num_args, char* args[]) {
	using namespace std;
	if (num_args <= 1) {
		cerr << "Need a problem name.\n";
		return 1;
	}
	string base_name = args[1];
	ifstream input(base_name + ".txt");
	const RetailProblem problem(input);
	problem.print_info();

	auto db = SqliteDb::fromFile("solutions.sqlite3");
	auto insert = db.make_statement<>(
	    "insert or ignore into solutions(name, objective, solution, solve_time, timestamp, "
	    "version, options)"
	    "VALUES(?1, ?2, ?3, ?4, strftime(\'%Y-%m-%dT%H:%M:%S\', \'now\'), ?5, ?6);");
	int best_value = numeric_limits<int>::max();
	double best_time = -1;
	auto callback = [&](const RetailLocalSearchInfo& info,
	                    const vector<vector<vector<int>>>& solution) {
		if (info.computed_objective_value < best_value) {
			best_value = info.computed_objective_value;
			best_time = info.elapsed_time;
		}
		insert.execute(base_name,
		               info.computed_objective_value,
		               problem.solution_to_string(solution),
		               info.elapsed_time,
		               version::revision,
		               "local_search " + get_all_command_line_options());
		return true;
	};
	retail_local_search(problem, callback);
	cout << base_name << " " << best_value << " " << best_time << endl;
	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
