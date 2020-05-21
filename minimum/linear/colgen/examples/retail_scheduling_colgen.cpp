#include <fstream>
#include <random>

#include <gflags/gflags.h>

#include <minimum/algorithms/dag.h>
#include <minimum/core/color.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/openmp.h>
#include <minimum/core/random.h>
#include <minimum/core/range.h>
#include <minimum/core/sqlite.h>
#include <minimum/linear/colgen/retail_scheduling_pricing.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/retail_scheduling.h>
#include <minimum/linear/solver.h>
#include <version/version.h>

using namespace std;
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

DEFINE_int32(num_solutions, 1, "Number of solutions to compute. Default: 1.");

class ShiftShedulingColgenProblem : public SetPartitioningProblem {
   public:
	ShiftShedulingColgenProblem(RetailProblem problem_)
	    : SetPartitioningProblem(problem_.staff.size(), problem_.num_cover_constraints()),
	      problem(move(problem_)) {
		for (int p : range(problem.staff.size())) {
			rng.emplace_back(repeatably_seeded_engine<std::mt19937>(p));
		}

		Timer timer("Setting up colgen problem");
		solution = make_grid<int>(problem.staff.size(), problem.periods.size(), problem.num_tasks);

		// Covering constraints with slack variables.
		for (auto d : range(problem.periods.size())) {
			auto& period = problem.periods[d];
			for (auto s : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(d, s);
				// Initialize this constraint.
				initialize_constraint(
				    c, period.min_cover.at(s), period.max_cover.at(s), 10000.0, 0, 0, 1.0);
			}
		}
		timer.OK();

		timer.next("Creating initial rosters");

		// Initial duals just try to cover as much as possible, with some
		// randomness to encourage different rosters.
		auto initial_duals = make_grid<double>(problem.staff.size(), number_of_rows());
		uniform_real_distribution<double> dual(-10, 100);
		for (int p = 0; p < problem.staff.size(); ++p) {
			for (int r = 0; r < number_of_rows(); ++r) {
				initial_duals[p][r] = dual(rng[p]);
			}

			// For testing, we can force staff to work a certain shift with the following code.
			int period = problem.day_start(0) + 2 * 4;
			int task = 0;
			int row = problem.staff.size() + problem.cover_constraint_index(period, task);
			initial_duals[p][row] = 1e6;
		}
		auto no_fixes =
		    make_grid<int>(problem.periods.size(), problem.num_tasks, []() { return -1; });

		OpenMpExceptionStore exception_store;
#pragma omp parallel for
		for (int p = 0; p < problem.staff.size(); ++p) {
			try {
				bool success = create_roster_graph(
				    problem, initial_duals[p], p, no_fixes, &solution[p], &rng[p]);
				check(success, "Initial roster generation failed.");
			} catch (...) {
				exception_store.store();
			}
		}
		exception_store.throw_if_available();

		for (int p = 0; p < problem.staff.size(); ++p) {
			add_column(create_column(p, solution[p]));
		}

		timer.OK();
	}

	~ShiftShedulingColgenProblem() {}

	Column create_column(int staff_index, const vector<vector<int>>& solution_for_staff) {
		const double cost = 0;

		Column column(cost, 0, 1);
		column.add_coefficient(staff_index, 1);

		const int start_row = problem.staff.size();
		for (auto d : range(problem.periods.size())) {
			for (auto s : range(problem.num_tasks)) {
				int c = problem.cover_constraint_index(d, s);
				int sol = solution_for_staff[d][s];
				minimum_core_assert(sol == 0 || sol == 1);
				if (sol == 1) {
					column.add_coefficient(start_row + c, 1.0);
				}
			}
		}
		return column;
	}

	virtual void generate(const std::vector<double>& dual_variables) override {
		// vector of ints instead of bool because we may be writing to
		// it concurrently.
		vector<int> success(problem.staff.size(), 0);

		OpenMpExceptionStore exception_store;
#pragma omp parallel for
		for (int p = 0; p < problem.staff.size(); ++p) {
			try {
				if (generate_for_staff(p, dual_variables)) {
					success[p] = 1;
				}
			} catch (...) {
				exception_store.store();
			}
		}
		exception_store.throw_if_available();

		for (int p = 0; p < problem.staff.size(); ++p) {
			if (success[p] == 1) {
				add_column(create_column(p, solution[p]));
			}
		}
	}

	virtual std::string member_name(int member) const override {
		return problem.staff.at(member).id;
	}

	std::vector<std::vector<std::vector<int>>> get_solution() {
		compute_fractional_solution(active_columns());

		for (int p = 0; p < problem.staff.size(); ++p) {
			for (int d = 0; d < problem.periods.size(); ++d) {
				for (int t = 0; t < problem.num_tasks; ++t) {
					int c = problem.cover_constraint_index(d, t);
					solution[p][d][t] = fractional_solution(p, c) > 0.5 ? 1 : 0;
				}
			}
		}

		return solution;
	}

   private:
	bool generate_for_staff(int staff_index, const std::vector<double>& dual_variables) {
		if (member_fully_fixed(staff_index)) {
			return false;
		}

		auto fixes = make_grid<int>(problem.periods.size(), problem.num_tasks);
		auto& fixes_for_staff = fixes_for_member(staff_index);
		minimum_core_assert(fixes_for_staff.size() == problem.num_cover_constraints());
		for (auto c : range(problem.num_cover_constraints())) {
			int d = problem.cover_constraint_to_period(c);
			int s = problem.cover_constraint_to_task(c);
			fixes[d][s] = fixes_for_staff[c];
		}

		return create_roster_graph(
		    problem, dual_variables, staff_index, fixes, &solution[staff_index], &rng[staff_index]);
	}

	const minimum::linear::RetailProblem problem;
	std::vector<std::vector<std::vector<int>>> solution;
	vector<mt19937> rng;
};

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
	if (num_args <= 1) {
		cerr << "Usage:\n";
		cerr << "    " << args[0] << " list\n ";
		cerr << "    " << args[0] << " run <base filename>\n ";
		cerr << "    " << args[0] << " print <base filename>\n ";
		cerr << "    " << args[0] << " save <base filename> <objective value>\n ";
		return 1;
	}
	string command = args[1];
	auto db = SqliteDb::fromFile("solutions.sqlite3");
	db.create_table("solutions",
	                {{"name", "string"},
	                 {"objective", "integer"},
	                 {"solution", "text"},
	                 {"solve_time", "real"},
	                 {"timestamp", "text"},
	                 {"version", "string"},
	                 {"options", "string"}});
	db.make_statement<>(
	      "create unique index if not exists objective_index on solutions (name, objective);")
	    .execute();

	if (command == "list") {
		auto select = db.make_statement<string, int>(
		    "select name, min(objective) from solutions "
		    "group by name "
		    "order by name;");
		for (auto& result : select.execute()) {
			cout << get<0>(result) << ": " << GREEN << get<1>(result) << NORMAL << ".\n";
		}
		return 0;
	}

	if (num_args <= 2) {
		cerr << "Need a base filename for " << command << ".\n";
		return 1;
	}

	string filename_base = args[2];
	ifstream input(filename_base + ".txt");
	const RetailProblem problem(input);
	problem.print_info();

	auto insert = db.make_statement<>(
	    "insert or ignore into solutions(name, objective, solution, solve_time, timestamp, "
	    "version, options)"
	    "VALUES(?1, ?2, ?3, ?4, strftime(\'%Y-%m-%dT%H:%M:%S\', \'now\'), ?5, ?6);");

	if (command == "print") {
		auto select = db.make_statement<int, string, double, string, string, string>(
		    "select * from ("
		    "  select objective, solution, solve_time, timestamp, version, ifnull(options, '') "
		    "  from solutions "
		    "  where name == ?1 "
		    "  order by objective asc "
		    "  limit 30 "
		    ") "
		    "order by objective desc;");
		for (auto& result : select.execute(filename_base)) {
			cout << YELLOW << setw(10) << right << get<0>(result) << NORMAL << ": "
			     << get<2>(result) << "s. Ran at " << get<3>(result) << " with (" << get<4>(result)
			     << "). " << get<5>(result) << "\n";
		}
	} else if (command == "save") {
		if (num_args <= 3) {
			cerr << "Need objetive value.\n";
			return 1;
		}
		int objective_value = from_string<int>(args[3]);
		auto select = db.make_statement<string, double, string>(
		    "select solution, solve_time, timestamp from solutions "
		    "where name == ?1 and objective == ?2;");
		auto result = select.execute(filename_base, objective_value).get();
		auto solution = problem.string_to_solution(get<0>(result));
		problem.save_solution(filename_base + ".ros",
		                      filename_base + ".xml",
		                      solution,
		                      get<1>(result),
		                      get<2>(result));
		cout << "Solution saved to " << filename_base << ".xml.\n";
	} else if (command == "run") {
		ShiftShedulingColgenProblem colgen_problem(problem);
		int best_objective_value = 1000'000'000;
		vector<vector<vector<int>>> best_solution;

		auto start_time = wall_time();
		for (int i = 1; i <= FLAGS_num_solutions; ++i) {
			colgen_problem.unfix_all();
			colgen_problem.solve();

			cerr << "Colgen done.\n";
			auto elapsed_time = wall_time() - start_time;
			cerr << "Elapsed time : " << elapsed_time << "s.\n";
			Timer timer("Saving solution...");
			auto solution = colgen_problem.get_solution();
			int objective_value = -1;
			try {
				objective_value = problem.check_feasibility(solution);
				if (objective_value < best_objective_value) {
					best_objective_value = objective_value;
					best_solution = solution;
				}
			} catch (std::runtime_error& error) {
				timer.fail();
				cerr << error.what() << endl;
			}
			if (objective_value >= 0) {
				insert.execute(filename_base,
				               objective_value,
				               problem.solution_to_string(solution),
				               elapsed_time,
				               version::revision,
				               get_all_command_line_options());
				timer.OK();
				cerr << "Objective value: " << objective_value << endl;
			}
		}

		if (!best_solution.empty()) {
			cerr << "Final objective value: " << best_objective_value << endl;
		}
	} else {
		cerr << "Unknown command " << command << "\n";
		return 2;
	}

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
