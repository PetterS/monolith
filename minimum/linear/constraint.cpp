
#include <vector>

#include <minimum/core/check.h>
#include <minimum/linear/constraint.h>
#include <minimum/linear/sum.h>

namespace minimum {
namespace linear {

class ConstraintList::Implementation {
   public:
	Implementation() {}

	template <typename T>
	Implementation(T&& t) : constraints(std::forward<T>(t)) {}

	std::vector<Constraint> constraints;
};

Constraint::Constraint(const LogicalExpression& expression)
    : lower_bound(1),
      upper_bound(static_cast<double>(expression.sum.size())),
      sum(expression.sum) {}

Constraint::Constraint(LogicalExpression&& expression)
    : lower_bound(1),
      upper_bound(static_cast<double>(expression.sum.size())),
      sum(std::move(expression.sum)) {}

Constraint::Constraint(double lower_bound_, Sum&& sum_, double upper_bound_)
    : lower_bound(lower_bound_), upper_bound(upper_bound_), sum(std::move(sum_)) {}

Constraint operator<=(const Variable& lhs, double rhs) { return {-1e100, lhs - rhs, 0}; }

Constraint operator<=(double lhs, const Variable& rhs) { return {-1e100, lhs - rhs, 0}; }

Constraint operator<=(const Variable& lhs, const Variable& rhs) { return {-1e100, lhs - rhs, 0}; }

Constraint operator<=(Sum lhs, const Sum& rhs) {
	lhs -= rhs;
	return Constraint(-1e100, std::move(lhs), 0.0);
}

Constraint operator<=(Sum lhs, double rhs) { return {-1e100, std::move(lhs), rhs}; }

Constraint operator<=(double lhs, Sum rhs) { return {lhs, std::move(rhs), 1e100}; }

Constraint operator>=(const Variable& lhs, double rhs) { return {0, lhs - rhs, 1e100}; }

Constraint operator>=(double lhs, const Variable& rhs) { return {0, lhs - rhs, 1e100}; }

Constraint operator>=(const Variable& lhs, const Variable& rhs) { return {0, lhs - rhs, 1e100}; }

Constraint operator>=(Sum lhs, const Sum& rhs) {
	lhs -= rhs;
	return Constraint(0.0, std::move(lhs), 1e100);
}

Constraint operator>=(Sum lhs, double rhs) { return {rhs, std::move(lhs), 1e100}; }

Constraint operator>=(double lhs, Sum rhs) { return {-1e100, std::move(rhs), lhs}; }

Constraint operator==(const Variable& lhs, double rhs) { return Constraint(0.0, lhs - rhs, 0.0); }

Constraint operator==(double lhs, const Variable& rhs) { return Constraint(0.0, lhs - rhs, 0.0); }

Constraint operator==(const Variable& lhs, const Variable& rhs) {
	return Constraint(0.0, lhs - rhs, 0.0);
}

Constraint operator==(Sum lhs, const Sum& rhs) {
	lhs -= rhs;
	return Constraint(0.0, std::move(lhs), 0.0);
}

Constraint operator==(Sum lhs, double rhs) { return {rhs, std::move(lhs), rhs}; }

Constraint operator==(double lhs, Sum rhs) { return {lhs, std::move(rhs), lhs}; }

ConstraintList::ConstraintList() : impl(new Implementation) {}

ConstraintList::ConstraintList(ConstraintList&& rhs) : impl(rhs.impl) { rhs.impl = nullptr; }

ConstraintList::ConstraintList(Constraint&& constraint) : impl(new Implementation) {
	impl->constraints.push_back(std::move(constraint));
}

ConstraintList& ConstraintList::operator&=(Constraint&& rhs) {
	impl->constraints.push_back(std::move(rhs));
	return *this;
}

ConstraintList::~ConstraintList() {
	if (impl) {
		delete impl;
	}
}

ConstraintList operator&&(Constraint&& lhs, Constraint&& rhs) {
	ConstraintList list(std::move(lhs));
	list &= std::move(rhs);
	return list;
}

ConstraintList operator&&(ConstraintList&& lhs, Constraint&& rhs) {
	ConstraintList list(std::move(lhs));
	list &= std::move(rhs);
	return list;
}

std::size_t ConstraintList::size() const { return impl->constraints.size(); }

const Constraint* ConstraintList::constraints() const {
	minimum_core_assert(size() > 0);
	return impl->constraints.data();
}
}  // namespace linear
}  // namespace minimum