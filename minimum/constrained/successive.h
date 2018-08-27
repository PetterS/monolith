#include <functional>
#include <optional>
#include <string_view>

#include <minimum/constrained/export.h>
#include <minimum/linear/glpk.h>
#include <minimum/nonlinear/constrained_function.h>

namespace minimum::constrained {
// Linearizes the objective and all constraints and solve the resulting
// linear program.
// Uses GLPK only because it can easily be warmstarted. If glp already
// contains a problem, presolve will be disabled and it will be updated.
// Otherwise, rows and columns will be added.
bool MINIMUM_CONSTRAINED_API
successive_linear_programming_step(minimum::linear::Solver* solver,
                                   std::optional<minimum::linear::SolutionsPointer>* solutions,
                                   minimum::nonlinear::ConstrainedFunction* constrained_function,
                                   double step_size);

class MINIMUM_CONSTRAINED_API SuccessiveLinearProgrammingSolver {
   public:
	std::function<void(std::string_view)> log_function = nullptr;
	double initial_step_size = 1;
	// Function improvement tolerance. The solver terminates
	// if the relative improvement is less than this and the
	// solution is feasible.
	double function_improvement_tolerance = 1e-3;
	int max_iterations = 300;

	enum class Result { OK = 0, NO_CONVERGENCE = 1, INFEASIBLE = 2 };
	Result solve(minimum::nonlinear::ConstrainedFunction* function);

   private:
};
}  // namespace minimum::constrained