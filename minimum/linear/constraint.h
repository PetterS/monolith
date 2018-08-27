#pragma once

#include <minimum/linear/sum.h>

#include <minimum/linear/export.h>

namespace minimum {
namespace linear {

// Represents a linear constraint.
//
// It is created from a sum and one of the
// operators <=, >= or ==.
//
//		auto c1 = x >= 0;
//		auto c2 = 3*x + 4*y == 4;
//
class MINIMUM_LINEAR_API Constraint {
	friend class IP;
	friend MINIMUM_LINEAR_API Constraint operator<=(const Variable& lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator<=(double lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator<=(const Variable& lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator<=(Sum lhs, const Sum& rhs);
	friend MINIMUM_LINEAR_API Constraint operator<=(Sum lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator<=(double lhs, Sum rhs);

	friend MINIMUM_LINEAR_API Constraint operator>=(const Variable& lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator>=(double lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator>=(const Variable& lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator>=(Sum lhs, const Sum& rhs);
	friend MINIMUM_LINEAR_API Constraint operator>=(Sum lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator>=(double lhs, Sum rhs);

	friend MINIMUM_LINEAR_API Constraint operator==(const Variable& lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator==(double lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator==(const Variable& lhs, const Variable& rhs);
	friend MINIMUM_LINEAR_API Constraint operator==(Sum lhs, const Sum& rhs);
	friend MINIMUM_LINEAR_API Constraint operator==(Sum lhs, double rhs);
	friend MINIMUM_LINEAR_API Constraint operator==(double lhs, Sum rhs);

   public:
	Constraint();  // For std::vector.
	Constraint(const LogicalExpression& expression);
	Constraint(LogicalExpression&& expression);

   private:
	Constraint(double lower_bound_, Sum&& sum_, double upper_bound_);
	double lower_bound, upper_bound;
	Sum sum;
};

MINIMUM_LINEAR_API Constraint operator<=(const Variable& lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator<=(double lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator<=(const Variable& lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator<=(Sum lhs, const Sum& rhs);
MINIMUM_LINEAR_API Constraint operator<=(Sum lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator<=(double lhs, Sum rhs);

MINIMUM_LINEAR_API Constraint operator>=(const Variable& lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator>=(double lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator>=(const Variable& lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator>=(Sum lhs, const Sum& rhs);
MINIMUM_LINEAR_API Constraint operator>=(Sum lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator>=(double lhs, Sum rhs);

MINIMUM_LINEAR_API Constraint operator==(const Variable& lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator==(double lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator==(const Variable& lhs, const Variable& rhs);
MINIMUM_LINEAR_API Constraint operator==(Sum lhs, const Sum& rhs);
MINIMUM_LINEAR_API Constraint operator==(Sum lhs, double rhs);
MINIMUM_LINEAR_API Constraint operator==(double lhs, Sum rhs);

class MINIMUM_LINEAR_API ConstraintList {
	friend class IP;

   public:
	ConstraintList(Constraint&& constraint);
	ConstraintList(ConstraintList&&);
	~ConstraintList();

	ConstraintList& operator&=(Constraint&&);

   private:
	ConstraintList();
	ConstraintList(const ConstraintList&);
	ConstraintList& operator=(const ConstraintList&) = delete;

	std::size_t size() const;
	const Constraint* constraints() const;

	class Implementation;
	Implementation* impl;
};

MINIMUM_LINEAR_API ConstraintList operator&&(Constraint&&, Constraint&&);
MINIMUM_LINEAR_API ConstraintList operator&&(ConstraintList&&, Constraint&&);
}  // namespace linear
}  // namespace minimum
