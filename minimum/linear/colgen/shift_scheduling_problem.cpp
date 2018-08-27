#include <fstream>
#include <random>
#include <vector>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/flamegraph.h>
#include <minimum/core/grid.h>
#include <minimum/core/openmp.h>
#include <minimum/core/range.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/shift_scheduling_pricing.h>
#include <minimum/linear/colgen/shift_scheduling_problem.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/scheduling_util.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;

DEFINE_string(pool_file_name,
              "",
              "If set, load column pool before starting and save from time to time. This is the "
              "prefix of the file name used. Default: not set.");

DEFINE_bool(save_solution, true, "Whether the solution should be saved to disk. Default: true.");

DEFINE_bool(best_integer_solution,
            false,
            "Use integer programming to compute solution. Very slow!");

namespace minimum {
namespace linear {
namespace colgen {

Column create_column(const minimum::linear::proto::SchedulingProblem& problem,
                     int p,
                     const vector<vector<int>>& solution_for_staff) {
	FLAMEGRAPH_LOG_FUNCTION;
	double cost = 0;

	// Shift preferences.
	for (const auto& pref : problem.worker(p).shift_preference()) {
		cost += pref.cost() * solution_for_staff[pref.day()][pref.shift()];
	}

	// Day-off preferences.
	for (const auto& pref : problem.worker(p).day_off_preference()) {
		bool working = false;
		for (auto s : range(problem.shift_size())) {
			if (solution_for_staff[pref.day()][s] == 1) {
				working = true;
			}
		}
		if (!working) {
			cost += pref.cost();
		}
	}

	Column column(cost, 0, 1);
	column.add_coefficient(p, 1);

	int r = problem.worker_size();
	for (const auto& requirement : problem.requirement()) {
		int sol = solution_for_staff[requirement.day()][requirement.shift()];
		minimum_core_assert(sol == 0 || sol == 1);
		if (sol == 1) {
			column.add_coefficient(r, 1.0);
		}
		++r;
	}
	return column;
}

struct ShiftShedulingColgenProblem::UserCommandData {
	// Graphs to save next iteration.  Will be cleared afterwards.
	// Maps staff ids to file names.
	map<string, string> save_graphs_next_iteration;
};

ShiftShedulingColgenProblem::ShiftShedulingColgenProblem(
    const minimum::linear::proto::SchedulingProblem& problem_)
    : SetPartitioningProblem(problem_.worker_size(), problem_.requirement_size()),
      problem(problem_),
      user_command_data(new UserCommandData) {
	FLAMEGRAPH_LOG_FUNCTION;

	// Append a unique number to the pool file name to be able to use the same prefix.
	if (FLAGS_pool_file_name != "") {
		pool_file_name = to_string(FLAGS_pool_file_name, ".", number_of_rows(), ".pool");
	}

	Timer timer("Setting up colgen problem");
	set_objective_constant(problem.objective_constant());

	// Covering constraints with slack variables.
	for (auto c : range(problem.requirement_size())) {
		initialize_constraint(c,
		                      problem.requirement(c).wanted(),
		                      problem.requirement(c).over_cost(),
		                      problem.requirement(c).under_cost());
	}

	day_shift_to_constraint =
	    make_grid<int>(problem.num_days(), problem.shift_size(), []() { return -1; });
	for (auto c : range(problem.requirement_size())) {
		auto& elem =
		    day_shift_to_constraint[problem.requirement(c).day()][problem.requirement(c).shift()];
		minimum_core_assert(elem == -1);
		elem = c;
	}
	// Verify that we have a cover constrait for all shifts.
	for (auto d : range(problem.num_days())) {
		for (auto s : range(problem.shift_size())) {
			minimum_core_assert(day_shift_to_constraint[d][s] >= 0);
		}
	}

	// Create scratch arrays.
	solution = make_grid<int>(problem.worker_size(), problem.num_days(), problem.shift_size());

	timer.OK();

	// Create and add initial bad roster to make problem feasible.
	timer.next("Creating initial rosters");

	// Initial duals just try to cover as much as possible, with some
	// randomness to encourage different rosters.
	auto initial_duals = make_grid<double>(problem.worker_size(), number_of_rows());
	uniform_real_distribution<double> dual(1, 100);
	for (int p = 0; p < problem.worker_size(); ++p) {
		for (int r = 0; r < number_of_rows(); ++r) {
			initial_duals[p][r] = dual(rng);
		}
	}

	auto no_fixes = make_grid<int>(problem.num_days(), problem.shift_size(), []() { return -1; });

	OpenMpExceptionStore exception_store;
#pragma omp parallel for
	for (int p = 0; p < problem.worker_size(); ++p) {
		try {
			minimum_core_assert(
			    create_roster_cspp(problem, initial_duals[p], p, no_fixes, &solution[p]));

		} catch (...) {
			exception_store.store();
		}
	}
	exception_store.throw_if_available();

	for (int p = 0; p < problem.worker_size(); ++p) {
		add_column(p, solution[p]);
	}
	timer.OK();

	cerr << "-- Objective for initial columns: ";
	try {
		cerr << objective_value(problem, solution);
	} catch (std::exception& e) {
		cerr << "Infeasible: " << e.what();
	}
	cerr << "\n";
}

ShiftShedulingColgenProblem::~ShiftShedulingColgenProblem() {}

bool ShiftShedulingColgenProblem::generate_for_staff(
    int p,
    const std::vector<double>& dual_variables,
    vector<vector<int>>* solution_for_staff) const {
	if (member_fully_fixed(p)) {
		return false;
	}

	auto fixes = make_grid<int>(problem.num_days(), problem.shift_size());
	auto& fixes_for_staff = fixes_for_member(p);
	for (auto c : range(problem.requirement_size())) {
		fixes[problem.requirement(c).day()][problem.requirement(c).shift()] = fixes_for_staff[c];
	}

	string graph_file_name;
	auto file_name_itr = user_command_data->save_graphs_next_iteration.find(problem.worker(p).id());
	if (file_name_itr != user_command_data->save_graphs_next_iteration.end()) {
		graph_file_name = file_name_itr->second;
	}

	return create_roster_cspp(
	    problem, dual_variables, p, fixes, solution_for_staff, graph_file_name);
}

double ShiftShedulingColgenProblem::integral_solution_value() {
	FLAMEGRAPH_LOG_FUNCTION;

	double objective;
	if (FLAGS_best_integer_solution) {
		auto ip = create_ip(true);
		GlpkSolver solver;
		solver.set_silent(true);
		minimum_core_assert(solver.solutions(ip.get())->get());
		objective = ip->get_entire_objective();
	} else {
		objective = SetPartitioningProblem::integral_solution_value();
	}

	return objective;
}

void ShiftShedulingColgenProblem::generate(const std::vector<double>& dual_variables) {
	FLAMEGRAPH_LOG_FUNCTION;
	// TOOD: Save column pool.
	// if (information.iteration % 100 == 0) {
	// 	possibly_save_column_pool();
	// }

	if (!loaded_pool_from_file) {
		if (pool_file_name != "") {
			try {
				Timer t(to_string("Loading pool from ", pool_file_name));
				load_column_pool();
				t.OK();
			} catch (std::runtime_error& e) {
				cerr << "-- Failed loading pool from " << pool_file_name << ":" << e.what() << "\n";
			}
		}
		loaded_pool_from_file = true;
	}

	// vector of ints instead of bool because we are writing to
	// it concurrently.
	vector<int> success(problem.worker_size(), 0);
	OpenMpExceptionStore exception_store;
#pragma omp parallel for
	for (int p = 0; p < problem.worker_size(); ++p) {
		try {
			if (generate_for_staff(p, dual_variables, &solution[p])) {
				success[p] = 1;
			}
		} catch (...) {
			exception_store.store();
		}
	}
	exception_store.throw_if_available();

	for (int p = 0; p < problem.worker_size(); ++p) {
		if (success[p] == 1) {
			// No need to waste a perfectly good column. If if is feasible for
			// other staff members, add it for them as well.
			for (int p2 = 0; p2 < problem.worker_size(); ++p2) {
				if (is_feasible_for_other(p, solution[p], p2)) {
					add_column(p2, solution[p]);
				}
			}
		}
	}

	user_command_data->save_graphs_next_iteration.clear();
}

const vector<vector<vector<int>>>& ShiftShedulingColgenProblem::get_solution() {
	compute_fractional_solution(active_columns());
	return get_solution_from_current_fractional();
}

const vector<vector<vector<int>>>&
ShiftShedulingColgenProblem::get_solution_from_current_fractional() {
	FLAMEGRAPH_LOG_FUNCTION;

	for (int p = 0; p < problem.worker_size(); ++p) {
		for (int d = 0; d < problem.num_days(); ++d) {
			for (int s = 0; s < problem.shift_size(); ++s) {
				solution[p][d][s] = fractional_solution(p, d, s) > 0.5 ? 1 : 0;
			}
		}
	}
	return solution;
}

void ShiftShedulingColgenProblem::possibly_save_column_pool() const {
	FLAMEGRAPH_LOG_FUNCTION;

	if (pool_file_name != "") {
		Timer t(to_string("Saving pool to ", pool_file_name));
		ofstream fout(pool_file_name, ios::binary);
		pool.save_to_stream(&fout);
		t.OK();
	}
}

void ShiftShedulingColgenProblem::load_column_pool() {
	if (pool_file_name != "") {
		ifstream fin(pool_file_name, ios::binary);
		ColumnPool new_pool(&fin);
		auto num_rows = number_of_rows();
		for (auto& column : new_pool) {
			int p = -1;
			for (auto& entry : column) {
				minimum_core_assert(0 <= entry.row && entry.row < num_rows, "Invalid pool file.");
				if (0 <= entry.row && entry.row < problem.worker_size()) {
					minimum_core_assert(entry.coef == 1, "Invalid pool file.");
					minimum_core_assert(p == -1, "Invalid pool file.");
					p = entry.row;
				}
			}
			if (p != -1) {
				SetPartitioningProblem::add_column(move(column));
			}
		}
	}
}

void ShiftShedulingColgenProblem::add_column(
    int p, const std::vector<std::vector<int>>& solution_for_staff) {
	SetPartitioningProblem::add_column(create_column(problem, p, solution_for_staff));
}

bool ShiftShedulingColgenProblem::interrupt_handler() {
	cerr << "\nCtrl-C\n";
	while (true) {
		cerr << ">";
		string command;
		cin.clear();
		getline(cin, command);

		auto args = split(command, ' ');
		if (args.empty()) {
		} else if (args.front() == "c" || args.front() == "continue") {
			return true;
		} else if (args.front() == "s" || args.front() == "stop") {
			throw std::runtime_error("User stop.");
		} else if (args.front() == "graph") {
			if (args.size() != 3) {
				cerr << "Invalid.\n";
			} else {
				user_command_data->save_graphs_next_iteration[args[1]] = args[2];
				cerr << "Will save graph for " << args[1] << " to " << args[2] << "\n";
			}
		} else {
			cerr << "(c)ontinue: continue optimization\n";
			cerr << "(s)stop: stops by throwing an exception\n";
			cerr << "graph <staff> <filename>: Saves the graph for this staff member next "
			        "iteration.\n";
		}
	}
	return false;
}

bool ShiftShedulingColgenProblem::is_feasible_for_other(int p,
                                                        const vector<vector<int>>& solution,
                                                        int p2) {
	if (p == p2) {
		return true;
	}
	if (problem.worker(p2).consecutive_shifts_limit().max()
	    < problem.worker(p).consecutive_shifts_limit().max()) {
		return false;
	}
	if (problem.worker(p2).consecutive_shifts_limit().min()
	    > problem.worker(p).consecutive_shifts_limit().min()) {
		return false;
	}
	if (problem.worker(p2).time_limit().max() < problem.worker(p).time_limit().max()) {
		return false;
	}
	if (problem.worker(p2).time_limit().min() > problem.worker(p).time_limit().min()) {
		return false;
	}
	if (problem.worker(p2).working_weekends_limit().max()
	    < problem.worker(p).working_weekends_limit().max()) {
		return false;
	}
	if (problem.worker(p2).working_weekends_limit().min()
	    > problem.worker(p).working_weekends_limit().min()) {
		return false;
	}
	if (problem.worker(p2).consecutive_days_off_limit().max()
	    < problem.worker(p).consecutive_days_off_limit().max()) {
		return false;
	}
	if (problem.worker(p2).consecutive_days_off_limit().min()
	    > problem.worker(p).consecutive_days_off_limit().min()) {
		return false;
	}
	for (int s : range(problem.shift_size())) {
		if (problem.worker(p2).shift_limit(s).max() < problem.worker(p).shift_limit(s).max()) {
			return false;
		}
		if (problem.worker(p2).shift_limit(s).min() > problem.worker(p).shift_limit(s).min()) {
			return false;
		}
	}
	for (int d : range(problem.num_days())) {
		for (int s : range(problem.shift_size())) {
			if (fix_state(p2, d, s) >= 0 && fix_state(p2, d, s) != solution[d][s]) {
				return false;
			}
		}
	}
	for (int d : problem.worker(p2).required_day_off()) {
		for (int s : range(problem.shift_size())) {
			if (solution[d][s] != 0) {
				return false;
			}
		}
	}
	return true;
}

int ShiftShedulingColgenProblem::fix_state(int p, int d, int s) const {
	int fix = fixes_for_member(p)[day_shift_to_constraint[d][s]];
	if (fix >= 0) {
		return fix;
	}
	for (int s2 : range(problem.shift_size())) {
		if (fixes_for_member(p)[day_shift_to_constraint[d][s2]] == 1) {
			return 0;
		}
	}
	return fix;
}

double ShiftShedulingColgenProblem::fractional_solution(int p, int d, int s) const {
	return SetPartitioningProblem::fractional_solution(p, day_shift_to_constraint[d][s]);
}
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
