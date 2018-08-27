#pragma once
#include <minimum/linear/solver.h>

class OsiSolverInterface;

namespace minimum {
namespace linear {

namespace constraint_solver {
class Searcher;
}

// New API for solving IP instances.
class MINIMUM_LINEAR_API ConstraintSolver : public Solver {
   public:
	ConstraintSolver() : Solver() {}
	~ConstraintSolver();

	virtual SolutionsPointer solutions(IP* ip) const override;
};
}  // namespace linear
}  // namespace minimum
