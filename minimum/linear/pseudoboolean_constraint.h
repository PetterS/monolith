#pragma once

#include <minimum/linear/variable.h>

#include <minimum/linear/export.h>
#include <minimum/linear/pseudoboolean.h>

namespace minimum {
namespace linear {

// Represents a pseudo-boolean constraint.
//
// It is created from a PseudoBoolean and one of the
// operators <=, >= or ==.
class MINIMUM_LINEAR_API PseudoBooleanConstraint {
	friend class IP;
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(PseudoBoolean lhs,
	                                                             PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(double lhs, PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(PseudoBoolean lhs, double rhs);

	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(PseudoBoolean lhs,
	                                                             PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(double lhs, PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(PseudoBoolean lhs, double rhs);

	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(PseudoBoolean lhs,
	                                                             PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(double lhs, PseudoBoolean rhs);
	friend MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(PseudoBoolean lhs, double rhs);

   private:
	PseudoBooleanConstraint(double lower_bound_, PseudoBoolean&& sum_, double upper_bound_);
	double lower_bound, upper_bound;
	PseudoBoolean sum;
};

MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(PseudoBoolean lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(double lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator<=(PseudoBoolean lhs, double rhs);

MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(PseudoBoolean lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(double lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator>=(PseudoBoolean lhs, double rhs);

MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(PseudoBoolean lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(double lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBooleanConstraint operator==(PseudoBoolean lhs, double rhs);
}  // namespace linear
}  // namespace minimum
