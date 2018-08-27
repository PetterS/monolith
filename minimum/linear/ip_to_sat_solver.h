#pragma once
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

class MINIMUM_LINEAR_API IpToSatSolver : public Solver {
   public:
	IpToSatSolver(std::function<std::unique_ptr<SatSolver>()> sat_solver_factory_);
	~IpToSatSolver();
	IpToSatSolver(const IpToSatSolver&) = delete;
	IpToSatSolver& operator=(const IpToSatSolver&) = delete;

	virtual SolutionsPointer solutions(IP* ip) const override;

	std::unique_ptr<SatSolver> convert(IP* ip);

	// The SAT solver can not use a cost function and will issue an
	// error if it is not 0. Calling this function allows the
	// solver to ignore the cost function if required and only
	// solve a feasibility problem.
	bool allow_ignoring_cost_function = false;

   private:
	std::function<std::unique_ptr<SatSolver>()> sat_solver_factory;
};
}  // namespace linear
}  // namespace minimum
