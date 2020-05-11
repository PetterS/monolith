#include <fstream>

#include <minimum/core/color.h>
#include <minimum/core/grid.h>
#include <minimum/core/main.h>
#include <minimum/core/range.h>
#include <minimum/linear/colgen/set_partitioning_problem.h>
#include <minimum/linear/retail_scheduling.h>

using namespace std;
using namespace minimum::core;
using namespace minimum::linear;
using namespace minimum::linear::colgen;

class ShiftShedulingColgenProblem : public SetPartitioningProblem {
   public:
	ShiftShedulingColgenProblem(RetailProblem problem_)
	    : SetPartitioningProblem(problem_.staff.size(), problem_.num_constraints()),
	      problem(move(problem_)) {
		Timer timer("Setting up colgen problem");
		solution = make_grid<int>(problem.staff.size(), problem.periods.size(), problem.num_tasks);

		// Covering constraints with slack variables.
		for (auto d : range(problem.periods.size())) {
			auto& period = problem.periods[d];
			for (auto s : range(problem.num_tasks)) {
				int c = constraint_index(d, s);
				// Initialize this constraint.
				initialize_constraint(
				    c, period.min_cover.at(s), period.max_cover.at(s), 1000.0, 1000.0);
			}
		}
		timer.OK();
	}

	~ShiftShedulingColgenProblem() {}

	virtual double integral_solution_value() override { return 0; }

	virtual void generate(const std::vector<double>& dual_variables) override {}

	virtual std::string member_name(int member) const override {
		return problem.staff.at(member).id;
	}

   private:
	int constraint_index(int period, int task) const { return period * problem.num_tasks + task; }

	const minimum::linear::RetailProblem problem;
	std::vector<std::vector<std::vector<int>>> solution;
};

int main_program(int num_args, char* args[]) {
	if (num_args <= 1) {
		cerr << "Need a file name.\n";
		return 1;
	}
	ifstream input(args[1]);
	const RetailProblem problem(input);
	problem.print_info();

	ShiftShedulingColgenProblem colgen_problem(problem);

	return 0;
}

int main(int num_args, char* args[]) { return main_runner(main_program, num_args, args); }
