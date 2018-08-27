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

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/sqlite.h>
using namespace minimum::core;

DEFINE_int32(limit, 5, "Limit of temporary and refined solutions to print. Default: 5.");

int main_program(int num_args, char* args[]) {
	using namespace std;
	check(num_args >= 2, "Need a file name.");
	auto db = SqliteDb::fromFile(args[1]);

	if (num_args >= 3) {
		int objective = from_string<int>(args[2]);
		for (auto table : {"solution", "refined", "temp"}) {
			auto select = db.make_statement<string_view>(
			    to_string("SELECT solution FROM ", table, " WHERE objective = ?1;"));
			auto result = select.execute(objective).possibly_get<0>();
			if (result.has_value()) {
				cout << *result << endl;
				return 0;
			}
		}
		std::cerr << "Solution not found.\n";
		return 1;
	} else {
		cout << "Here are the best temporary solutions from colgen:\n";
		auto select_temp = db.make_statement<int, string_view>(
		    "SELECT * FROM "
		    "  (SELECT objective, options "
		    "  FROM temp "
		    "  ORDER BY objective ASC "
		    "  LIMIT ?1) "
		    "ORDER BY objective DESC;");
		for (auto& result : select_temp.execute(FLAGS_limit)) {
			cout << YELLOW << setw(7) << right << get<0>(result) << NORMAL << " with "
			     << get<1>(result) << "\n";
		}

		try {
			auto select_refined = db.make_statement<int, int, double>(
			    "SELECT * FROM "
			    "  (SELECT objective, start_objective, time "
			    "  FROM refined "
			    "  ORDER BY objective ASC "
			    "  LIMIT ?1) "
			    "ORDER BY objective DESC;");
			cout << "Here are the best refined solutions:\n";
			for (auto& result : select_refined.execute(FLAGS_limit)) {
				cout << YELLOW << setw(7) << right << get<0>(result) << NORMAL << " from original "
				     << get<1>(result) << " in " << get<2>(result) << "s.\n";
			}
		} catch (std::runtime_error&) {
			// Table does not have to exist.
		}

		cout << "Here are all colgen solutions:\n";
		auto select_solution = db.make_statement<int, int, double, string_view>(
		    "SELECT objective, number, time, options "
		    "FROM solution "
		    "ORDER BY objective DESC, number ASC;");
		for (auto& result : select_solution.execute()) {
			cout << YELLOW << setw(7) << right << get<0>(result) << NORMAL << " in " << WHITE
			     << get<2>(result) << "s" << NORMAL << " as number " << get<1>(result) << " with "
			     << get<3>(result) << "\n";
		}
	}
	return 0;
}

int main(int num_args, char* args[]) {
	if (num_args >= 2 && ends_with(to_string(args[1]), "-help")) {
		gflags::SetUsageMessage(
		    "Runs colgen method for shift scheduling.\n   shift_scheduling_colgen [options] "
		    "<input_file>");
		gflags::ShowUsageWithFlagsRestrict("shift_scheduling_colgen", "colgen");
		return 0;
	}
	gflags::FlagSaver saver;  // Because Emscripten calls main multiple times.
	gflags::ParseCommandLineFlags(&num_args, &args, true);
	try {
		return main_program(num_args, args);
	} catch (std::exception& err) {
		cerr << "Error: " << err.what() << endl;
		return 1;
	}
}
