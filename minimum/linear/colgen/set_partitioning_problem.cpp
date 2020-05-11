#include <cmath>
#include <vector>
using namespace std;

#include <gflags/gflags.h>

#include <minimum/core/check.h>
#include <minimum/core/flamegraph.h>
#include <minimum/core/grid.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/proto.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
using namespace minimum::core;

DEFINE_double(fix_threshold,
              0.75,
              "Treshold to use when fixing constraints to members. Default: 0.75.");

DEFINE_double(
    objective_change_before_fixing,
    1.0,
    "If the objective function decreases by more than this, fixing will not run. Default: 1.0.");

DEFINE_double(zero_fix_threshold,
              0.05,
              "Treshold to use when fixing constraints to zero for a member. Default: 0.05.");

namespace minimum {
namespace linear {
namespace colgen {

class SetPartitioningProblem::Implementation {
   public:
	Implementation(int number_of_groups_,
	               int number_of_constraints_,
	               SetPartitioningProblem* parent_)
	    : number_of_groups(number_of_groups_),
	      number_of_constraints(number_of_constraints_),
	      constraints(number_of_constraints_),
	      columns_for_member(number_of_groups_),
	      num_fixes_for_member(number_of_groups_, 0),
	      parent(parent_) {
		fixes = make_grid<int>(number_of_groups, number_of_constraints, []() { return -1; });
		fractional_solution = make_grid<double>(number_of_groups, number_of_constraints);
		// Start with a history of 0.5.
		fractional_solution_history =
		    make_grid<double>(number_of_groups, number_of_constraints, []() { return 0.5; });
	}

	// Check that invariants hold.
	void sanity_check() const {
		for (int p : range(number_of_groups)) {
			for (auto i : columns_for_member[p]) {
				minimum_core_assert(column_member[i] == p);
			}
		}

		for (auto i : range(parent->pool.size())) {
			int p = column_member[i];
			if (p < 0) {
				// Not a member column.
				continue;
			}
			int found = 0;
			for (auto i2 : columns_for_member[p]) {
				if (i == i2) {
					found++;
				}
			}
			minimum_core_assert(found++);
		}
	}

	// Associates column i with member p.
	void set_column_member(size_t i, int p) {
		if (column_member.size() <= i) {
			column_member.resize(i + 1, -1);
		}
		column_member[i] = p;
		columns_for_member.at(p).push_back(i);
	}

	void print_fractionality_stats() const {
		cerr << "-- Members remaining: ";
		int total_remaining = 0;
		for (int p = 0; p < number_of_groups; ++p) {
			int nonzeros = 0;
			for (auto i : columns_for_member[p]) {
				auto& column = parent->pool.at(i);
				if (column.solution_value > 1e-6) {
					nonzeros++;
				}
			}
			minimum_core_assert(nonzeros > 0);
			if (nonzeros > 1) {
				cerr << "(" << parent->member_name(p) << ", " << nonzeros << ") ";
				total_remaining++;
			}
		}
		cerr << "and total is " << total_remaining << ".\n";
	}

	void update_fractional_history() {
		constexpr double update_fraction = 0.95;
		for (int p = 0; p < number_of_groups; ++p) {
			for (int c : range(number_of_constraints)) {
				fractional_solution_history[p][c] =
				    update_fraction * fractional_solution_history[p][c]
				    + (1.0 - update_fraction) * fractional_solution[p][c];
			}
		}
	}

	// Calls fix_constraint_to_member repeatedly for all members. Returns the
	// number of columns fixes.
	int fix_constraints(double threshold, const std::vector<std::size_t>& active_columns) {
		FLAMEGRAPH_LOG_FUNCTION;

		int fixed = 0;
		vector<string> fixed_for_members;
		vector<string> fixed_for_members_to_zero;
		for (int p = 0; p < number_of_groups; ++p) {
			int member_fixes = 0;
			int member_fixes_to_zero = 0;
			for (auto c : range(number_of_constraints)) {
				if (fixes[p][c] < 0 && fractional_solution[p][c] >= threshold) {
					fixed += fix_constraint_to_member(p, c, 1);

					// Consider the following three columns
					// 0: fractional 0.4,  A  .  C
					// 1: fractional 0.2,  A  B  .
					// 2: fractional 0.4,  .  B  C
					//
					// If we first fix C, 1 is removed. If we then fix B, 0 is
					// removed. We must therefore recompute the fractional solution
					// after each fix.
					parent->compute_fractional_solution(p, active_columns);

					member_fixes++;
				} else if (fixes[p][c] < 0
				           && fractional_solution_history[p][c] <= FLAGS_zero_fix_threshold) {
					fix_constraint_to_member(p, c, 0);
					// Just as above, recompute the fractional solution.
					parent->compute_fractional_solution(p, active_columns);
					member_fixes_to_zero++;
				}
			}
			if (member_fixes > 0) {
				fixed_for_members.push_back(to_string(parent->member_name(p), ": ", member_fixes));
			}
			if (member_fixes_to_zero > 0) {
				fixed_for_members_to_zero.push_back(
				    to_string(parent->member_name(p), ": ", member_fixes_to_zero));
			}
		}
		if (!fixed_for_members.empty()) {
			cerr << "-- Fixed for " << to_string(fixed_for_members) << " at level " << threshold
			     << ".\n";
			proto::Event event;
			event.mutable_fix_shift()->set_level(threshold);
			parent->attach_event(move(event));
		}
		if (!fixed_for_members_to_zero.empty()) {
			cerr << "-- Fixed to zero for " << to_string(fixed_for_members_to_zero) << " at level "
			     << FLAGS_zero_fix_threshold << ".\n";
		}

		// Only return the number fixed to one.
		return fixed;
	}

	// Adds a fix that a group member always contributes 0 or 1 to a specific
	// constraint.
	// After calling this function with value = 1, the member must have
	// that constraint in all of their new columns.
	// After calling this function with value = 0, the member can not have
	// that constraint in any of their new columns.
	//
	// Returns the number of columns fixed (to zero, fixing always results in
	// columns not fulfilling the criteria becoming zero).
	int fix_constraint_to_member(int member, int constraint, int value) {
		minimum_core_assert(value == 0 || value == 1);
		if (fixes[member][constraint] == value) {
			return 0;
		}
		minimum_core_assert(fixes[member][constraint] < 0);
		fixes[member][constraint] = value;
		num_fixes_for_member[member]++;

		int fixed_columns = 0;
		bool at_least_one_column_left = false;
		for (auto i : columns_for_member[member]) {
			auto& column = parent->pool.at(i);
			minimum_core_assert(column_member[i] == member);

			if (column.is_fixed()) {
				if (column.upper_bound() > 0.5) {
					at_least_one_column_left = true;
				}
				continue;
			}

			bool column_has_constraint = false;
			for (auto& entry : column) {
				if (entry.row >= number_of_groups) {
					int c = entry.row - number_of_groups;
					if (c == constraint) {
						column_has_constraint = true;
						break;
					}
				}
			}

			if (value == 1) {
				if (column_has_constraint) {
					at_least_one_column_left = true;
				} else {
					fixed_columns++;
					column.fix(0);
				}
			} else {
				if (column_has_constraint) {
					fixed_columns++;
					column.fix(0);
				} else {
					at_least_one_column_left = true;
				}
			}
		}
		minimum_core_assert(at_least_one_column_left, "Problem will become infeasible.");

		return fixed_columns;
	}

	// Rounds the columns for all groups and stores the integer solution to be retrieved by
	// fractional_solution(p, i). Does not modify the solution value of the columns.
	//
	// The rounding is done by picking the column with the highest fractional value for
	// each group. Thus, the solution will be feasible if picking whole columns always
	// results in a feasible solution. Rounding in a more clever way (as done when fixing)
	// could give an infeasible solution.
	//
	// Returns the cost of all columns for group members.
	double make_solution_integer(const std::vector<std::size_t>& active_columns) {
		vector<pair<size_t, double>> best_column(number_of_groups, make_pair(0, 0));
		vector<pair<size_t, double>> saved_solution_values;

		for (auto i : active_columns) {
			auto p = column_member[i];
			if (p < 0) {
				continue;
			}
			auto& column = parent->pool.at(i);
			if (column.solution_value > best_column[p].second) {
				best_column[p].first = i;
				best_column[p].second = column.solution_value;
			}

			saved_solution_values.emplace_back(i, column.solution_value);
			column.solution_value = 0;
		}

		double column_cost = 0;
		for (auto& entry : best_column) {
			minimum_core_assert(entry.second > 0);
			column_cost += parent->pool.at(entry.first).cost();
			parent->pool.at(entry.first).solution_value = 1;
		}

		// Compute and store the solution from the columns we just rounded.
		parent->compute_fractional_solution(active_columns);

		// Restore the fractional solution.
		for (auto& entry : saved_solution_values) {
			parent->pool.at(entry.first).solution_value = entry.second;
		}

		return column_cost;
	}

	// Computes the cost from the constraints from the solution stored in
	// fractional_solution.
	double cover_cost(const std::vector<std::size_t>& active_columns) const {
		vector<double> rows(number_of_constraints, 0);
		for (auto member : range(number_of_groups)) {
			for (auto c : range(number_of_constraints)) {
				rows[c] += fractional_solution[member][c];
			}
		}

		double cost = 0;
		for (auto c : range(number_of_constraints)) {
			if (rows[c] > constraints[c].max_value) {
				cost += constraints[c].over_cost * (rows[c] - constraints[c].max_value);
			} else if (rows[c] < constraints[c].min_value) {
				cost += constraints[c].under_cost * (constraints[c].min_value - rows[c]);
			}
		}
		return cost;
	}

	int number_of_groups;
	int number_of_constraints;

	struct Constraint {
		double min_value = numeric_limits<double>::quiet_NaN();
		double max_value = numeric_limits<double>::quiet_NaN();
		double over_cost = numeric_limits<double>::quiet_NaN();
		double under_cost = numeric_limits<double>::quiet_NaN();
	};
	vector<Constraint> constraints;

	// The member for this column. Can be -1 for e.g. slack columns.
	vector<int> column_member;
	// All columns for a given member.
	vector<vector<size_t>> columns_for_member;

	vector<vector<int>> fixes;
	vector<int> num_fixes_for_member;
	vector<vector<double>> fractional_solution;
	vector<vector<double>> fractional_solution_history;

	SetPartitioningProblem* parent;
};

SetPartitioningProblem::SetPartitioningProblem(int number_of_groups_, int number_of_constraints_)
    : Problem(number_of_groups_ + number_of_constraints_),
      impl(new Implementation(number_of_groups_, number_of_constraints_, this)) {
	for (auto m : range(number_of_groups_)) {
		set_row_lower_bound(m, 1);
		set_row_upper_bound(m, 1);
	}
}

SetPartitioningProblem::~SetPartitioningProblem() { delete impl; }

std::string SetPartitioningProblem::member_name(int member) const { return to_string(member); }

int SetPartitioningProblem::column_member(int column) const {
	return impl->column_member.at(column);
}

int SetPartitioningProblem::number_of_rows() const {
	int num_rows = 0;
	// One SOS1 constraint per group.
	num_rows += impl->number_of_groups;
	// The rest of the constraints come after.
	num_rows += impl->number_of_constraints;
	return num_rows;
}

void SetPartitioningProblem::add_column(Column&& column) {
	check(column.lower_bound() == 0 && column.upper_bound() == 1,
	      "SetPartitioningProblem::add_column: Only 0/1 columns accepted.");

	int p = -1;
	for (auto& entry : column) {
		if (entry.row < impl->number_of_groups) {
			check(p < 0, "Invalid column.");
			check(entry.coef == 1.0, "Invalid column.");
			p = entry.row;
		} else {
			check(entry.coef == 1.0,
			      "SetPartitioningProblem::add_column: Column needs to have a 1 as coefficient for "
			      "all "
			      "constraints. The fixing logic currently assumes that.");
		}
	}
	check(p >= 0,
	      "Columns with no member (e.g. slack columns) should be added to the pool directly.");

	auto s1 = pool.size();
	pool.add(std::move(column));
	auto s2 = pool.size();
	if (s1 < s2) {
		impl->set_column_member(s1, p);
	}
}

int SetPartitioningProblem::fix(const FixInformation& information) {
	FLAMEGRAPH_LOG_FUNCTION;
	return fix_using_columns(information, active_columns());
}

int SetPartitioningProblem::fix_using_columns(const FixInformation& information,
                                              const std::vector<std::size_t>& active_columns) {
	FLAMEGRAPH_LOG_FUNCTION;

	impl->sanity_check();
	compute_fractional_solution(active_columns);
	impl->update_fractional_history();

	if (information.iteration <= 4
	    || information.objective_change < -FLAGS_objective_change_before_fixing
	    || information.objective_change > 0) {
		return 0;
	}

	int fixed = 0;
	double threshold = FLAGS_fix_threshold;
	// Iteratively lower the threshold if we aren't able to fix
	// any constraints.
	while (fixed == 0 && threshold > 0.51) {
		fixed = impl->fix_constraints(threshold, active_columns);
		threshold *= 0.9;
	}

	// If no constraints were fixed, fix a column to a member. This is a much
	// stronger fix, as this member will not be optimized any more after this.
	if (fixed == 0) {
		for (auto i : active_columns) {
			auto& column = pool.at(i);
			if (column.upper_bound() == 1 && column.solution_value >= 0.5
			    && column.solution_value < 1.0 - 1e-7) {
				int p = impl->column_member[i];

				cerr << "-- Fixing column " << i << " to 1 for member " << member_name(p)
				     << ", was " << column.solution_value << "\n";
				proto::Event event;
				event.mutable_fix_column()->set_fractional_value(column.solution_value);
				attach_event(move(event));

				fixed++;
				impl->num_fixes_for_member[p] = impl->number_of_constraints;
				column.fix(1);
				break;
			}
		}
	}

	// At this point, just fix anything to proceed.
	if (fixed == 0) {
		for (int p = 0; p < impl->number_of_groups && fixed == 0; ++p) {
			int nonzeros = 0;
			for (auto i : impl->columns_for_member[p]) {
				auto& column = pool.at(i);
				if (column.solution_value > 1e-6) {
					nonzeros++;
				}
			}
			if (nonzeros > 1) {
				for (auto i : impl->columns_for_member[p]) {
					auto& column = pool.at(i);
					if (!column.is_fixed() && column.solution_value > 1e-6) {
						cerr << "-- Fixing column " << i << " to 1 for member " << member_name(p)
						     << " as a last resort, was " << column.solution_value << "\n";
						proto::Event event;
						event.mutable_fix_column()->set_fractional_value(column.solution_value);
						attach_event(move(event));

						column.fix(1);
						fixed++;
						impl->num_fixes_for_member[p] = impl->number_of_constraints;
						for (auto j : impl->columns_for_member[p]) {
							if (i != j) {
								pool.at(j).fix(0);
								fixed++;
							}
						}
						break;
					}
				}
			}
		}
	}

	impl->print_fractionality_stats();

	if (fixed == 0) {
		// Just try to round the existing fractional solution to an integer one.
		cerr << "-- Not able to fix anything. Creating integer solution.\n";
		// This solution is hopefully feasible. We are done.
		// fixes = solution;
		fixed = -1;
	}

	return fixed;
}

void SetPartitioningProblem::initialize_constraint(
    int c, double min_value, double max_value, double over_cost, double under_cost) {
	impl->constraints.at(c).min_value = min_value;
	impl->constraints.at(c).max_value = max_value;
	impl->constraints.at(c).over_cost = over_cost;
	impl->constraints.at(c).under_cost = under_cost;

	int r = impl->number_of_groups + c;
	set_row_lower_bound(r, min_value);
	set_row_upper_bound(r, max_value);

	Column under_slack(under_cost, 0, 1e100);
	under_slack.add_coefficient(r, 1.0);
	under_slack.set_integer(false);
	pool.add(move(under_slack));

	Column over_slack(over_cost, 0, 1e100);
	over_slack.add_coefficient(r, -1.0);
	over_slack.set_integer(false);
	pool.add(move(over_slack));
}

void SetPartitioningProblem::unfix_all() {
	for (int p = 0; p < impl->number_of_groups; ++p) {
		for (int i = 0; i < impl->number_of_constraints; ++i) {
			impl->fixes[p][i] = -1;
			impl->fractional_solution_history[p][i] = 0.5;
		}
	}
	for (auto& column : pool) {
		column.unfix();
	}
	fill(impl->num_fixes_for_member.begin(), impl->num_fixes_for_member.end(), 0);
}

bool SetPartitioningProblem::member_fully_fixed(int member) const {
	return impl->num_fixes_for_member.at(member) >= impl->number_of_constraints;
}

double SetPartitioningProblem::integral_solution_value(
    const std::vector<std::size_t>& active_columns) {
	auto column_cost = impl->make_solution_integer(active_columns);
	return column_cost + impl->cover_cost(active_columns);
}

double SetPartitioningProblem::integral_solution_value() {
	return integral_solution_value(active_columns());
}

void SetPartitioningProblem::compute_fractional_solution(
    const std::vector<std::size_t>& active_columns) {
	FLAMEGRAPH_LOG_FUNCTION;
	for (int p = 0; p < impl->number_of_groups; ++p) {
		compute_fractional_solution(p, active_columns);
		for (auto c : range(impl->number_of_constraints)) {
			if (fixes_for_member(p)[c] >= 0) {
				minimum_core_assert(
				    abs(fixes_for_member(p)[c] - fractional_solution(p, c)) <= 1e-3,
				    "Fix state does not correspond to fractional solution. A column "
				    "violating the fix state may have been added.");
			}
		}
	}
}

void SetPartitioningProblem::compute_fractional_solution(
    int member, const std::vector<std::size_t>& active_columns) {
	for (auto c : range(impl->number_of_constraints)) {
		impl->fractional_solution[member][c] = 0;
	}

	// Compute the sum over all active columns for this member
	// in order to force the sum to equal 1.
	double member_sum = 0.0;
	for (auto i : active_columns) {
		auto& column = pool.at(i);
		if (impl->column_member[i] != member) {
			continue;
		}
		member_sum += column.solution_value;
	}

	for (auto i : active_columns) {
		auto& column = pool.at(i);
		if (impl->column_member[i] != member) {
			continue;
		}

		// The member sum should be exactly 1.0. Two things may cause it to deviate:
		// an incorrect first-order solver and the fact that this member function
		// is also called after performing a fix.
		minimum_core_assert(member_sum > 1e-6, "Member sum is ", member_sum);
		column.solution_value /= member_sum;
		column.solution_value = max(column.lower_bound(), column.solution_value);
		column.solution_value = min(column.upper_bound(), column.solution_value);

		auto eps = 1e-6;
		minimum_core_assert(column.lower_bound() - eps <= column.solution_value
		                        && column.solution_value <= column.upper_bound() + eps,
		                    column.lower_bound(),
		                    " <= ",
		                    column.solution_value,
		                    " <= ",
		                    column.upper_bound());

		for (auto& entry : column) {
			if (entry.row >= impl->number_of_groups) {
				int c = entry.row - impl->number_of_groups;
				impl->fractional_solution[member][c] += column.solution_value;
			}
		}
	}
}

bool SetPartitioningProblem::column_allowed(int column) const {
	int member = column_member(column);
	if (member < 0) {
		return true;
	}
	auto& fixes = fixes_for_member(member);
	vector<int> found(fixes.size(), 0);
	for (auto& entry : pool.at(column)) {
		if (entry.row >= impl->number_of_groups) {
			int c = entry.row - impl->number_of_groups;
			if (fixes.at(c) == 0) {
				return false;
			}
			found[c] = 1;
		}
	}
	for (auto i : range(fixes.size())) {
		if (fixes[i] == 1 && found[i] != 1) {
			return false;
		}
	}
	return true;
}

double SetPartitioningProblem::fractional_solution(int member, int constraint) const {
	return impl->fractional_solution.at(member).at(constraint);
}

const std::vector<int>& SetPartitioningProblem::fixes_for_member(int member) const {
	return impl->fixes.at(member);
}

proto::LogEntry SetPartitioningProblem::create_log_entry() const {
	proto::LogEntry entry;
	auto set_partitioning_entry = entry.mutable_set_partitioning();
	for (int p = 0; p < impl->number_of_groups; ++p) {
		auto member = set_partitioning_entry->add_member();
		member->set_id(member_name(p));

		for (auto c : range(impl->number_of_constraints)) {
			auto constraint = member->add_constraint();
			if (fixes_for_member(p)[c] == 0) {
				constraint->set_fixed(true);
				constraint->set_value(0);
			} else if (fixes_for_member(p)[c] == 1) {
				constraint->set_fixed(true);
				constraint->set_value(1);
			} else {
				constraint->set_value(impl->fractional_solution[p][c]);
			}
		}
	}
	return entry;
}
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
