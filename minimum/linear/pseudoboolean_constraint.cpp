
#include <minimum/core/check.h>
#include <minimum/linear/pseudoboolean_constraint.h>
using minimum::core::check;

namespace minimum {
namespace linear {

PseudoBooleanConstraint::PseudoBooleanConstraint(double lower_bound_,
                                                 PseudoBoolean&& sum_,
                                                 double upper_bound_)
    : lower_bound(lower_bound_), upper_bound(upper_bound_), sum(std::move(sum_)) {}

PseudoBooleanConstraint operator<=(PseudoBoolean lhs, PseudoBoolean rhs) {
	return {-1e100, std::move(lhs) - std::move(rhs), 0};
}

PseudoBooleanConstraint operator<=(double lhs, PseudoBoolean rhs) {
	return {lhs, std::move(rhs), 1e100};
}

PseudoBooleanConstraint operator<=(PseudoBoolean lhs, double rhs) {
	return {-1e100, std::move(lhs), rhs};
}

PseudoBooleanConstraint operator>=(PseudoBoolean lhs, PseudoBoolean rhs) {
	return {0, std::move(lhs) - std::move(rhs), 1e100};
}

PseudoBooleanConstraint operator>=(double lhs, PseudoBoolean rhs) {
	return {-1e100, std::move(rhs), lhs};
}

PseudoBooleanConstraint operator>=(PseudoBoolean lhs, double rhs) {
	return {rhs, std::move(lhs), 1e100};
}

PseudoBooleanConstraint operator==(PseudoBoolean lhs, PseudoBoolean rhs) {
	return {0, std::move(lhs) - std::move(rhs), 0};
}

PseudoBooleanConstraint operator==(double lhs, PseudoBoolean rhs) {
	return {lhs, std::move(rhs), lhs};
}

PseudoBooleanConstraint operator==(PseudoBoolean lhs, double rhs) {
	// std::clog << lhs << " == " << rhs << "\n";
	return {rhs, std::move(lhs), rhs};
}
}  // namespace linear
}  // namespace minimum
