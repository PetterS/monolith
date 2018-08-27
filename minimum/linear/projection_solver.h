#pragma once

#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {
// A very, very simple solver that just projects onto the
// constraints repeatedly.
// https://web.stanford.edu/class/ee392o/alt_proj.pdf
class MINIMUM_LINEAR_API ProjectionSolver : public Solver {
   public:
	virtual SolutionsPointer solutions(IP* ip) const override;
	double tolerance = 1e-4;
	int max_iterations = 100'000;
	int test_interval = 10'000;

   private:
};
}  // namespace linear
}  // namespace minimum
