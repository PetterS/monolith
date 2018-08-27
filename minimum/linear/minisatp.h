#pragma once
#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

class MINIMUM_LINEAR_API MinisatpSolver : public Solver {
   public:
	virtual SolutionsPointer solutions(IP* ip) const override;
};
}  // namespace linear
}  // namespace minimum
