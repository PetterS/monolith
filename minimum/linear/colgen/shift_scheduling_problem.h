#pragma once
#include <random>
#include <vector>

#include <gflags/gflags.h>

#include <minimum/linear/colgen/column.h>
#include <minimum/linear/colgen/export.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/proto.h>

DECLARE_string(pool_file_name);
DECLARE_bool(save_solution);

namespace minimum {
namespace linear {
namespace colgen {

MINIMUM_LINEAR_COLGEN_API minimum::linear::colgen::Column create_column(
    const minimum::linear::proto::SchedulingProblem& problem,
    int p,
    const std::vector<std::vector<int>>& solution_for_staff);

class MINIMUM_LINEAR_COLGEN_API ShiftShedulingColgenProblem : public SetPartitioningProblem {
   public:
	ShiftShedulingColgenProblem(const minimum::linear::proto::SchedulingProblem& problem_);
	~ShiftShedulingColgenProblem();

	// Generate a column for a single staff member. This member function
	// is const to ensure thread-safety.
	bool generate_for_staff(int p,
	                        const std::vector<double>& dual_variables,
	                        std::vector<std::vector<int>>* solution_for_staff) const;

	virtual double integral_solution_value() override;

	virtual void generate(const std::vector<double>& dual_variables) override;

	// Fills in solution with the current rounded solution from the
	// active columns.
	const std::vector<std::vector<std::vector<int>>>& get_solution();

	const std::vector<std::vector<std::vector<int>>>& get_solution_from_current_fractional();

	void possibly_save_column_pool() const;

	void load_column_pool();

	virtual std::string member_name(int member) const override {
		return problem.worker(member).id();
	}

	void add_column(int p, const std::vector<std::vector<int>>& solution_for_staff);

   protected:
	virtual bool interrupt_handler() override;

	const minimum::linear::proto::SchedulingProblem& problem;

   private:
	// Whether a feasible solution for p is also feasible for p2.
	bool is_feasible_for_other(int p, const std::vector<std::vector<int>>& solution, int p2);

	int fix_state(int p, int d, int s) const;
	double fractional_solution(int p, int d, int s) const;

	std::string pool_file_name;

	std::vector<std::vector<std::vector<int>>> solution;

	mutable std::vector<std::mt19937_64> random_engines;

	std::vector<std::vector<int>> day_shift_to_constraint;

	bool loaded_pool_from_file = false;

	std::mt19937_64 rng;

	struct UserCommandData;
	// Used to save commands entered by the user after interruption.
	std::unique_ptr<UserCommandData> user_command_data;
};
}  // namespace colgen
}  // namespace linear
}  // namespace minimum
