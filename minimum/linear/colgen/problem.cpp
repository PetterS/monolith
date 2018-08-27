
#include <cmath>
#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <set>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/flamegraph.h>
#include <minimum/core/numeric.h>
#include <minimum/core/range.h>
#include <minimum/core/record_stream.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/time.h>
#include <minimum/linear/colgen/problem.h>
#include <minimum/linear/first_order_solver.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scs.h>
#include <minimum/linear/solver.h>
using namespace minimum::core;

DEFINE_int32(min_rmp_iterations,
             4,
             "Minimum number of iterations to run before fixing or stopping");

DEFINE_bool(use_first_order_solver,
            false,
            "Use a first-order solver instead of simplex. Default: false.");

DEFINE_bool(disable_rescale_c, false, "Rescales c in Chambolle-Pock. Default: false.");

DEFINE_bool(use_scs, false, "Uses SCS instead of Chambolle-Pock. Default: false.");

DEFINE_string(save_ip_directory,
              "",
              "When set, specifies the path prefix where all linear programes will be saved. "
              "Default: not set.");

DEFINE_string(
    proto_log_file,
    "",
    "Logs extra iteration data, e.g. the fractional solution for each member when fixing.");

namespace {
volatile std::sig_atomic_t signal_raised = 0;
void set_signal_raised(int signal) {
	signal_raised = 1;
	// Strange, but this is required in Visual Studio 2015 in order
	// to make signals work.
	// http://stackoverflow.com/questions/41552081/sigint-handler-reset-in-visual-c-2015
	std::signal(SIGINT, set_signal_raised);
}
}  // namespace

namespace minimum {
namespace linear {
namespace colgen {

constexpr int column_inactive_limit = 5;

class Problem::Implementation {
   public:
	Implementation(Problem& parent_, int number_of_rows_)
	    : parent(parent_),
	      number_of_rows(number_of_rows_),
	      row_lb(number_of_rows_, -1e100),
	      row_ub(number_of_rows_, 1e100),
	      dual_solution(number_of_rows, 0) {
		if (!FLAGS_proto_log_file.empty()) {
			iteration_data.emplace(FLAGS_proto_log_file, ios::binary);
		}
	}

	Problem& parent;
	const int number_of_rows;
	size_t initial_column_count = 0;
	vector<double> row_lb;
	vector<double> row_ub;
	vector<std::size_t> active_columns;
	vector<std::size_t> new_active_columns;
	vector<double> dual_solution;
	vector<int> column_inactive_count;
	double objective_constant = 0;

	std::optional<std::ofstream> iteration_data;
	std::vector<proto::Event> current_events;

	// Pick columns to include in the next linear program.
	void pick_columns(int iteration) {
		FLAMEGRAPH_LOG_FUNCTION;
		if (active_columns.empty()) {
			column_inactive_count.clear();
			for (auto i : range(initial_column_count)) {
				active_columns.push_back(i);
			}
			return;
		}
		if (iteration == 1) {
			column_inactive_count.clear();
			// First iteration, but we have active columns. This is
			// a subsequent run and we do not need to price in new
			// columns in the first round.
			return;
		}

		column_inactive_count.resize(parent.pool.size(), 0);

		new_active_columns.clear();
		// Any existing columns that are non-zero needs to still
		// be present in order to guarantee feasibility. Drop all
		// columns that are zero.
		for (auto i : active_columns) {
			auto& column = parent.pool.at(i);
			minimum_core_assert(column.solution_value >= -1e-5);

			if (column.solution_value < 1e-6) {
				// Before fixing, and with an exact solver, this column
				// will not have a negative reduced cost.
				column_inactive_count[i]++;
			} else {
				column_inactive_count[i] = 0;
			}

			// Keep a column if it was one of the initial columns or it has been active recently
			// and is not fixed.
			bool initial_column = i < initial_column_count;
			bool recently_active = column_inactive_count[i] < column_inactive_limit;
			bool fixed_to_zero = column.is_fixed() && column.upper_bound() == 0;
			if (initial_column || (recently_active && !fixed_to_zero)) {
				new_active_columns.push_back(i);
			}
		}

		// Place a (large) limit on the number of new columns. Otherwise there will be
		// an enormous amount added after restarting.
		auto new_good_columns = parent.pool.get_sorted(dual_solution, active_columns.size());
		for (auto& entry : new_good_columns) {
			new_active_columns.push_back(entry.index);
		}

		sort(new_active_columns.begin(), new_active_columns.end());
		auto end = unique(new_active_columns.begin(), new_active_columns.end());
		new_active_columns.erase(end, new_active_columns.end());

		active_columns.swap(new_active_columns);
	}

	std::unique_ptr<IP> create_ip(bool use_integer_variables,
	                              vector<Variable>* column_variables,
	                              vector<DualVariable>* row_variables) const {
		FLAMEGRAPH_LOG_FUNCTION;

		auto ip = std::make_unique<IP>();
		vector<Sum> row_sums(number_of_rows);
		vector<bool> row_used(number_of_rows, false);
		column_variables->clear();
		row_variables->clear();

		for (auto i : active_columns) {
			const auto& column = parent.pool.at(i);

			if (use_integer_variables && column.is_integer()) {
				column_variables->push_back(ip->add_variable(IP::Integer));
			} else {
				column_variables->push_back(ip->add_variable(IP::Real));
			}

			check(column.lower_bound() >= 0, "Columns can not be allowed to be negative.");
			ip->add_bounds(column.lower_bound(), column_variables->back(), column.upper_bound());
			ip->add_objective(column.cost() * column_variables->back());

			for (auto& entry : column) {
				row_sums.at(entry.row) += entry.coef * column_variables->back();
				row_used[entry.row] = true;
			}

			ip->set_solution(column_variables->back().get_index(), column.solution_value);
		}

		auto zero_variable_hack =
		    ip->add_variable(IP::Real);  // To avoid non-existing dual variables.
		ip->add_bounds(0, zero_variable_hack, 0);
		ip->set_solution(zero_variable_hack.get_index(), 0);

		for (int i = 0; i < number_of_rows; ++i) {
			minimum_core_assert(row_used[i]);
			row_variables->push_back(
			    ip->add_constraint(row_lb[i], row_sums[i] + zero_variable_hack, row_ub[i]));

			ip->set_dual_solution(row_variables->back().get_index(), dual_solution[i]);
		}

		if (!FLAGS_save_ip_directory.empty()) {
			Timer t("Saving IP");
			static int ip_count = 0;
			auto file_name = to_string(FLAGS_save_ip_directory, "/", ip_count++, ".ip");
			ofstream fout(file_name, ios::binary);
			check(ip->get().SerializeToOstream(&fout), "Save IP failed.");
			t.OK();
		}

		return ip;
	}

	double solve_lp(bool use_integer_variables = false) {
		FLAMEGRAPH_LOG_FUNCTION;
		vector<Variable> column_variables;
		vector<DualVariable> row_variables;
		auto ip = create_ip(use_integer_variables, &column_variables, &row_variables);

		check(!FLAGS_use_scs || FLAGS_use_first_order_solver,
		      "--scs requires --use_first_order_solver.");
		if (FLAGS_use_first_order_solver) {
			if (FLAGS_use_scs) {
				ScsSolver solver;
				solver.set_convergence_tolerance(1e-3);
				solver.set_silent(true);
				solver.set_max_iterations(5'000);

				solver.solutions(ip.get())->get();
			} else {
				PrimalDualSolver solver;
				solver.options().maximum_iterations = 5000;
				solver.options().tolerance = 1e-6;
				solver.options().rescale_c = !FLAGS_disable_rescale_c;

				solver.solutions(ip.get())->get();
			}
		} else {
			IPSolver solver;
			// GlopSolver solver;

			solver.set_silent(true);
			if (!solver.solutions(ip.get())->get()) {
				return numeric_limits<double>::quiet_NaN();
			}
		}

		double objective = ip->get_entire_objective();

		for (int i = 0; i < number_of_rows; ++i) {
			dual_solution[i] = row_variables[i].value();
		}

		for (auto j : range(active_columns.size())) {
			auto i = active_columns.at(j);
			auto lb = parent.pool.at(i).lower_bound();
			auto ub = parent.pool.at(i).upper_bound();
			parent.pool.at(i).solution_value = min(ub, max(lb, column_variables[j].value()));
		}

		return objective;
	}

	void write_proto_log_entry(proto::LogEntry log_entry) {
		at_scope_exit(current_events.clear());
		if (!iteration_data) {
			return;
		}
		for (auto& event : current_events) {
			auto new_event = log_entry.add_event();
			*new_event = move(event);
		}
		write_record(&*iteration_data, log_entry);
	}
};

Problem::Problem(int number_of_rows)
    : proto_logging_enabled(!FLAGS_proto_log_file.empty()),
      impl(new Implementation(*this, number_of_rows)) {}

Problem::~Problem() { delete impl; }

double Problem::integral_solution_value() { return numeric_limits<double>::quiet_NaN(); }

proto::LogEntry Problem::create_log_entry() const { return {}; }

void Problem::attach_event(proto::Event event) { impl->current_events.emplace_back(move(event)); }

double Problem::solve() {
	FLAMEGRAPH_LOG_FUNCTION;

	auto old_handler = std::signal(SIGINT, set_signal_raised);
	at_scope_exit(if (old_handler != SIG_ERR) { std::signal(SIGINT, old_handler); });
	signal_raised = 0;

	double global_start_time = wall_time();

	double objective = numeric_limits<double>::quiet_NaN();
	double previous_objective = numeric_limits<double>::quiet_NaN();
	double objective_change = numeric_limits<double>::quiet_NaN();

	if (impl->initial_column_count == 0) {
		impl->initial_column_count = pool.size();
	}

	for (int iteration = 1; iteration <= 100'000; ++iteration) {
		int generated_columns = 0;
		int fixed_columns = 0;

		double fix_time = numeric_limits<double>::quiet_NaN();
		double generate_time = numeric_limits<double>::quiet_NaN();
		double solve_time = numeric_limits<double>::quiet_NaN();

		if (iteration > FLAGS_min_rmp_iterations) {
			double start_time = wall_time();
			FixInformation information;
			information.iteration = iteration;
			information.objective_change = objective_change;
			fixed_columns = fix(information);
			fix_time = wall_time() - start_time;
		}

		if (iteration >= 2) {
			double start_time = wall_time();
			auto old_size = pool.size();
			generate(impl->dual_solution);
			generated_columns = pool.size() - old_size;
			generate_time = wall_time() - start_time;
		}

		double start_time = wall_time();
		auto previous_active_size = impl->active_columns.size();
		impl->pick_columns(iteration);
		auto active_size_change =
		    ptrdiff_t(impl->active_columns.size()) - ptrdiff_t(previous_active_size);

		objective = impl->solve_lp() + impl->objective_constant;
		objective_change = objective - previous_objective;
		previous_objective = objective;
		solve_time = wall_time() - start_time;

		size_t fractional_integer_active_columns = 0;
		size_t integer_active_columns = 0;
		for (auto i : impl->active_columns) {
			auto& column = pool.at(i);
			if (column.is_integer()) {
				integer_active_columns++;
				if (integer_residual(column.solution_value) > 1e-5) {
					fractional_integer_active_columns++;
				}
			}
		}

		proto::LogEntry log_entry;
		if (impl->iteration_data) {
			log_entry = create_log_entry();
		}

		if (iteration % 100 == 1) {
			// clang-format off
			cerr << "  Iter |       Objective      |    Problem   |    Pool   | Fixed  |                  Time                  | Frac. |\n";
			cerr << "       |    frac.     rounded |   size  diff | size gen. |        |    fix       gen     solve   tot. cum. |       |\n";
			// clang-format on
		}

		start_time = wall_time();
		auto rounded_value = integral_solution_value() + impl->objective_constant;
		string rounded_value_string = "n/a";
		if (rounded_value == rounded_value) {
			rounded_value_string = to_string(setprecision(0), fixed, rounded_value);
		}
		fix_time += wall_time() - start_time;

		log_entry.set_iteration(iteration);
		log_entry.set_fractional_objective(objective);
		log_entry.set_integer_objective(rounded_value);
		log_entry.set_active_size(impl->active_columns.size() - impl->initial_column_count);
		log_entry.set_active_size_change(active_size_change);
		log_entry.set_pool_size(pool.allowed_size());
		log_entry.set_pool_size_change(generated_columns);
		log_entry.set_fixed_columns(fixed_columns);
		log_entry.set_cumulative_time(wall_time() - global_start_time);

		cerr << to_string(setw(7), right, iteration) << " "
		     << to_string(setw(10), right, log_entry.fractional_objective()) << " "
		     << to_string(setw(10), right, rounded_value_string) << " "
		     << to_string(setw(8), right, log_entry.active_size()) << " "
		     << to_string(setw(4), right, log_entry.active_size_change()) << " "
		     << to_string(setw(7), right, log_entry.pool_size()) << " "
		     << to_string(setw(4), right, log_entry.pool_size_change()) << " "
		     << to_string(setw(6), right, log_entry.fixed_columns()) << " "
		     << to_string(setw(11), right, setprecision(1), scientific, fix_time) << " "
		     << to_string(setw(9), right, setprecision(2), scientific, generate_time) << " "
		     << to_string(setw(9), right, setprecision(2), scientific, solve_time) << " "
		     << to_string(setw(10), right, setprecision(3), scientific, log_entry.cumulative_time())
		     << " "
		     << to_string(setw(6),
		                  right,
		                  setprecision(0),
		                  fixed,
		                  100.0 * double(fractional_integer_active_columns)
		                      / double(integer_active_columns),
		                  "%")
		     << "\n";
		const bool integer_solution = iteration > FLAGS_min_rmp_iterations
		                              && abs(objective_change) <= 1e-6
		                              && fractional_integer_active_columns == 0;

		if (integer_solution) {
			proto::Event event;
			event.mutable_integer_solution();
			minimum_core_assert(event.has_integer_solution());
			attach_event(move(event));
		}
		impl->write_proto_log_entry(move(log_entry));

		if (objective != objective) {
			check(false, "Solving main LP failed.");
		}

		if (integer_solution) {
			cerr << "-- Stopping. Objective change is " << objective_change
			     << " with an integer feasible solution.\n";
			break;
		}
		if (fixed_columns < 0) {
			cerr << "-- Stopping. Requested by subclass.\n";
			break;
		}

		if (signal_raised) {
			if (!interrupt_handler()) {
				break;
			}
			signal_raised = 0;
		}
	}
	return objective;
}

const std::vector<size_t>& Problem::active_columns() const { return impl->active_columns; }

void Problem::set_objective_constant(double constant) { impl->objective_constant = constant; }

bool Problem::interrupt_handler() {
	throw std::runtime_error("Ctrl-C caught.");
	return false;
}

std::unique_ptr<IP> Problem::create_ip(bool use_integer_variables) const {
	vector<Variable> column_variables;
	vector<DualVariable> row_variables;
	return impl->create_ip(use_integer_variables, &column_variables, &row_variables);
}

void Problem::set_row_lower_bound(int row, double lower) { impl->row_lb.at(row) = lower; }

void Problem::set_row_upper_bound(int row, double upper) { impl->row_ub.at(row) = upper; }
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
