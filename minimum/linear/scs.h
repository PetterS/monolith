#pragma once
#include <string>

#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

class MINIMUM_LINEAR_API ScsSolver : public Solver {
   public:
	ScsSolver() : Solver() {}
	virtual ~ScsSolver();

	virtual SolutionsPointer solutions(IP* ip) const override;

	void set_max_iterations(int max);

	void set_convergence_tolerance(double eps);

   private:
	int max_iterations = -1;
	double convergence_tolerance = -1;
};
}  // namespace linear
}  // namespace minimum
