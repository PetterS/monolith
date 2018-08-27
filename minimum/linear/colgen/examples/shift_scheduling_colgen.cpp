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
#include <minimum/core/concurrent.h>
#include <minimum/core/flamegraph.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/openmp.h>
#include <minimum/core/range.h>
#include <minimum/core/sqlite.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/problem.h>
#include <minimum/linear/colgen/shift_scheduling_pricing.h>
#include <minimum/linear/colgen/shift_scheduling_problem.h>
#include <minimum/linear/data/util.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/scheduling_util.h>
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

DEFINE_int32(num_solutions, 1, "Number of solutions to compute. Default: 1.");

DEFINE_bool(benchmark_pricing, false, "Runs a small pricing benchmark and exits. Default: false.");

DEFINE_int32(verify_solution_value,
             -1,
             "If non-negative, verifies that all solutions obtained are better than this value. "
             "Default: -1.");

DEFINE_bool(read_gent_instance,
            false,
            "Instead loads an instance from "
            "http://www.projectmanagement.ugent.be/?q=research/personnel_scheduling/nsp");

DEFINE_bool(load_solutions,
            false,
            "If set, load columns from previous solutions in DB. Default: false.");

class MyShiftShedulingColgenProblem : public ShiftShedulingColgenProblem {
   public:
	MyShiftShedulingColgenProblem(const minimum::linear::proto::SchedulingProblem& problem_,
	                              std::string file_name_)
	    : ShiftShedulingColgenProblem(problem_), db(SqliteDb::inMemory()) {
		FLAMEGRAPH_LOG_FUNCTION;

		if (FLAGS_save_solution) {
			auto output_file_name = to_string(file_name_, ".db");
			db = SqliteDb::fromFile(output_file_name);
		}
		db.make_statement(
		      "CREATE TABLE IF NOT EXISTS "
		      "solution ("
		      "objective INTEGER UNIQUE,"
		      "solution TEXT,"
		      "number INTEGER,"
		      "time FLOAT,"
		      "options TEXT"
		      ");")
		    .execute();
		db.make_statement(
		      "CREATE TABLE IF NOT EXISTS "
		      "temp ("
		      "objective INTEGER UNIQUE,"
		      "solution TEXT,"
		      "options TEXT"
		      ");")
		    .execute();
		insert_solution =
		    db.make_statement("INSERT OR IGNORE INTO solution VALUES(?1, ?2, ?3, ?4, ?5);");
		insert_temp = db.make_statement("INSERT OR IGNORE INTO temp VALUES(?1, ?2, ?3);");
		smallest_temp = db.make_statement<int>("SELECT MIN(objective) FROM temp;");
	}

	double integral_solution_value() override {
		FLAMEGRAPH_LOG_FUNCTION;

		auto objective = ShiftShedulingColgenProblem::integral_solution_value();
		if (objective < best_integer_solution) {
			best_integer_solution = objective;
			auto solution = get_solution_from_current_fractional();
			write_solution(solution, -1, -1);
		}
		return objective;
	}

	void write_solution(const vector<vector<vector<int>>>& solution,
	                    int solution_index,
	                    double elapsed_time) {
		FLAMEGRAPH_LOG_FUNCTION;

		// Do not write all temp solutions.
		if (solution_index < 0) {
			auto smallest = smallest_temp->execute().get<0>();
			if (smallest != 0 && best_integer_solution + problem.objective_constant() >= smallest) {
				return;
			}
		}

		stringstream sout_solution;
		stringstream sout_flags;
		stringstream null;
		auto objective = print_solution(null, problem, solution);
		save_solution(sout_solution, problem, solution);
		vector<gflags::CommandLineFlagInfo> all_flags;
		gflags::GetAllFlags(&all_flags);
		for (auto option : all_flags) {
			if (!option.is_default) {
				sout_flags << option.name << "=" << option.current_value << " ";
			}
		}

		if (solution_index < 0) {
			insert_temp->execute(objective, sout_solution.str(), sout_flags.str());
		} else {
			insert_solution->execute(
			    objective, sout_solution.str(), solution_index, elapsed_time, sout_flags.str());
		}
	}

	void load_previous_solutions() {
		Timer t("Loading previous solution from DB");
		for (auto table : {"solution", "temp"}) {
			auto data_stmt =
			    db.make_statement<string>(to_string("SELECT solution FROM ", table, ";"));
			for (auto data : data_stmt.execute()) {
				istringstream solution_file(get<0>(data));
				auto loaded_solution = load_solution(solution_file, problem);
				for (int p = 0; p < problem.worker_size(); ++p) {
					add_column(p, loaded_solution[p]);
				}
			}
		}
		t.OK();
	}

   private:
	SqliteDb db;
	optional<SqliteStatement<>> insert_solution;
	optional<SqliteStatement<>> insert_temp;
	optional<SqliteStatement<int>> smallest_temp;

	double best_integer_solution = std::numeric_limits<double>::max();
};

int main_program(int num_args, char* args[]) {
	using namespace std;
	flamegraph::enable();
	check(num_args >= 2, "Need a file name.");

	auto start_time = wall_time();
	minimum::linear::proto::SchedulingProblem problem;
	if (FLAGS_read_gent_instance) {
		auto dir = data::get_directory();
		ifstream problem_file(dir + "/NSPLib/N25/1.nsp");
		ifstream case_file(dir + "/NSPLib/Cases/6.gen");
		problem = read_gent_instance(problem_file, case_file);
	} else {
		ifstream fin(args[1]);
		problem = read_nottingham_instance(fin);
		fin.close();
	}
	verify_scheduling_problem(problem);

	if (FLAGS_benchmark_pricing) {
		mt19937_64 rng;
		vector<double> initial_duals(problem.worker_size() + problem.requirement_size(), 100.0);
		uniform_real_distribution<double> dual_dist(1, 100);
		for (auto& dual : initial_duals) {
			dual = dual_dist(rng);
		}

		auto minus_one = []() { return -1; };
		auto fixes = make_grid<int>(problem.num_days(), problem.shift_size(), minus_one);
		auto solution = make_grid<int>(problem.num_days(), problem.shift_size());
		Timer t("Benchmarking pricing for one staff member");
		minimum_core_assert(create_roster_cspp(problem, initial_duals, 0, fixes, &solution));
		t.OK();
		return 0;
	}

	MyShiftShedulingColgenProblem colgen_problem(problem, args[1]);
	if (FLAGS_load_solutions) {
		colgen_problem.load_previous_solutions();
	}

	vector<double> solution_objectives;
	for (int i = 1; i <= FLAGS_num_solutions; ++i) {
		colgen_problem.unfix_all();
		auto problem_objective = colgen_problem.solve();
		auto& solution = colgen_problem.get_solution();
		cerr << "Colgen done.\n";
		auto elapsed_time = wall_time() - start_time;

		colgen_problem.possibly_save_column_pool();

		auto objective = print_solution(cout, problem, solution);
		minimum_core_assert(abs(problem_objective - objective) / (abs(objective) + 1e-5) <= 1e-4,
		                    "Objective from colgen does not match the computed one.");

		cerr << "Elapsed time : " << elapsed_time << "s.\n";
		solution_objectives.push_back(objective);

		if (FLAGS_verify_solution_value >= 0) {
			check(objective <= FLAGS_verify_solution_value,
			      "Solution of ",
			      objective,
			      " is not better than ",
			      FLAGS_verify_solution_value);
		}

		colgen_problem.write_solution(solution, i, elapsed_time);
	}

	if (FLAGS_num_solutions > 1) {
		cerr << "Solution objectives: " << to_string(solution_objectives) << ".\n";
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
	return main_runner(main_program, num_args, args);
}
