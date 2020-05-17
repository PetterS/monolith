#pragma once
#include <vector>

#include <gflags/gflags.h>

#include <minimum/linear/colgen/problem.h>

DECLARE_double(fix_threshold);
DECLARE_double(objective_change_before_fixing);

namespace minimum {
namespace linear {
namespace colgen {

class MINIMUM_LINEAR_COLGEN_API SetPartitioningProblem : public Problem {
   public:
	SetPartitioningProblem(int number_of_groups_, int number_of_constraints_);
	SetPartitioningProblem(const SetPartitioningProblem&) = delete;
	~SetPartitioningProblem();

	virtual int fix(const FixInformation& information) override;

	// Initializes a constraint. To be called once for every constraint before
	// solving starts.
	//
	// If the row of the constraint is more than requirement. The cost is
	//    (row - requirement) * over_cost.
	// Similarily for under_cost.
	//
	// This function will add two slack columns to the pool, which will then
	// be treated as inital columns.
	void initialize_constraint(int c,
	                           double min_value,
	                           double max_value,
	                           double under_cost,
	                           double over_cost,
	                           double under_quadratic_cost = 0,
	                           double over_quadratic_cost = 0);

	// Removes all perfomed fixes so the the problem can be solved again
	// for a new solution.
	void unfix_all();

	// Whether this member is fully determined by fixes. No need to
	// generate any more columns.
	bool member_fully_fixed(int member) const;

	// Rounds the current solution to an integer solution and returns its
	// cost.
	virtual double integral_solution_value() override;

	// Override to display a name for a member instead of an index.
	virtual std::string member_name(int member) const;

	// The member associated with a specific column. The member is the first
	// non-zero row of the column and this method is just a convenience.
	//
	// Note that not all columns correspond to members and this function
	// returns -1 in that case.
	int column_member(int column) const;

	virtual proto::LogEntry create_log_entry() const override;

   protected:
	int number_of_rows() const;

	// Visible for testing. could instead refactor fix() to take the active
	// columns.
	int fix_using_columns(const FixInformation& information,
	                      const std::vector<std::size_t>& active_columns);

	// Adds a column for a group member to the pool. This function should
	// be called instead of adding to the pool directly.
	//
	// Other columns (e.g. slack columns) should be added to the pool
	// directly on construction.
	//
	// TODO: This class should take care of the set covering constraints
	//       for the columns.
	void add_column(Column&& column);

	// Returns all fixes for a particular member. The returned vector has the
	// same size as the number of constraints. A negative value means no fix
	// for that particular constraint.
	const std::vector<int>& fixes_for_member(int member) const;

	// Whether a column in the pool is valid with respect to the fix state.
	bool column_allowed(int c) const;

	// How much a particular member contributes to a certain constraint.
	// Available after computing fractional solution.
	double fractional_solution(int member, int constraint) const;
	void compute_fractional_solution(const std::vector<std::size_t>& active_columns);
	void compute_fractional_solution(int member, const std::vector<std::size_t>& active_columns);

	// Version of integral_solution_value that takes the active columns as an argument, which
	// makes testing easier.
	double integral_solution_value(const std::vector<std::size_t>& active_columns);

   private:
	class Implementation;
	Implementation* impl;
};
}  // namespace colgen
}  // namespace linear
}  // namespace minimum