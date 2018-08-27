#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/main.h>
#include <minimum/core/record_stream.h>
#include <minimum/linear/colgen/proto.h>
using namespace minimum::core;

int main_program(int num_args, char* args[]) {
	using namespace std;
	check(num_args >= 2, "Need a file name.");

	ifstream file(args[1], ios::binary);
	check(bool{file}, "Could not open file.");

	while (file.peek() != EOF) {
		auto log_entry = read_record<minimum::linear::colgen::proto::LogEntry>(&file);
		int iteration = log_entry.iteration();
		if (iteration % 100 == 1) {
			// clang-format off
			cout << "  Iter |       Objective      |    Problem   |    Pool   | Fixed  |  Time      |\n";
			cout << "       |    frac.     rounded |   size  diff | size gen. |        |  tot. cum. |\n";
			// clang-format on
		}
		string rounded_value_string = "n/a";
		if (log_entry.integer_objective() == log_entry.integer_objective()) {
			rounded_value_string = to_string(setprecision(0), fixed, log_entry.integer_objective());
		}
		cout << to_string(setw(7), right, iteration) << " "
		     << to_string(setw(10), right, log_entry.fractional_objective()) << " "
		     << to_string(setw(10), right, rounded_value_string) << " "
		     << to_string(setw(8), right, log_entry.active_size()) << " "
		     << to_string(setw(4), right, log_entry.active_size_change()) << " "
		     << to_string(setw(7), right, log_entry.pool_size()) << " "
		     << to_string(setw(4), right, log_entry.pool_size_change()) << " "
		     << to_string(setw(8), right, log_entry.fixed_columns()) << "  "
		     << to_string(setw(10), right, setprecision(3), scientific, log_entry.cumulative_time())
		     << "\n";
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
