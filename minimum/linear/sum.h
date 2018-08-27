#pragma once

#include <minimum/linear/export.h>
#include <minimum/linear/variable.h>

namespace minimum {
namespace linear {

// A weighted sum of several variables (plus a constant).
//
// It is created by the user in several ways:
//
//		Sum  s1 = 0;
//		Sum  s2 = x;
//		auto s3 = x + y;
//		auto s4 = 2.0*x - y + 3;
//
class MINIMUM_LINEAR_API Sum {
	friend class Constraint;
	friend class IP;
	friend class LogicalExpression;
	friend class PseudoBoolean;

   public:
	Sum();
	Sum(const Sum& sum);
	Sum& operator=(const Sum& sum);
	Sum(Sum&& sum);
	Sum& operator=(Sum&& sum);
	Sum(double constant_);
	Sum(const Variable& variable);
	~Sum();

	Sum& operator+=(const Sum& rhs);
	Sum& operator-=(const Sum& rhs);
	Sum& operator*=(double coeff);
	Sum& operator/=(double coeff);
	void negate();

	// The value of the sum after the corresponding IP
	// has been solved.
	double value() const;

   protected:
	class Implementation;
	Implementation* impl;

	// The number of terms in the sum (not including constants).
	std::size_t size() const;
	// Indices of the variables in the sum (size() must be > 0
	// in order to call this method.
	const int* indices() const;
	// Coefficients of the variables in the sum (size() > 0).
	const double* values() const;
	// Constant term in the sum.
	double constant() const;
	// The IP that created the sum (or nullptr).
	const IP* creator() const;

	void match_solvers(const Sum& sum);
};

MINIMUM_LINEAR_API Sum operator*(double coeff, const Variable& variable);
MINIMUM_LINEAR_API Sum operator*(const Variable& variable, double coeff);
MINIMUM_LINEAR_API Sum operator*(double coeff, Sum sum);
MINIMUM_LINEAR_API Sum operator*(Sum sum, double coeff);

MINIMUM_LINEAR_API Sum operator/(Sum sum, double coeff);

MINIMUM_LINEAR_API Sum operator+(Variable lhs, Variable rhs);
MINIMUM_LINEAR_API Sum operator+(Variable lhs, Sum rhs);
MINIMUM_LINEAR_API Sum operator+(Sum lhs, Variable rhs);
MINIMUM_LINEAR_API Sum operator+(const Sum& lhs, const Sum& rhs);
MINIMUM_LINEAR_API Sum operator+(const Sum& lhs, Sum&& rhs);
MINIMUM_LINEAR_API Sum operator+(Sum&& lhs, const Sum& rhs);
MINIMUM_LINEAR_API Sum operator+(Sum&& lhs, Sum&& rhs);

MINIMUM_LINEAR_API Sum operator-(const Variable& variable);
MINIMUM_LINEAR_API Sum operator-(Sum rhs);

MINIMUM_LINEAR_API Sum operator-(Variable lhs, Variable rhs);
MINIMUM_LINEAR_API Sum operator-(Variable lhs, Sum rhs);
MINIMUM_LINEAR_API Sum operator-(Sum lhs, Variable rhs);
MINIMUM_LINEAR_API Sum operator-(const Sum& lhs, const Sum& rhs);
MINIMUM_LINEAR_API Sum operator-(const Sum& lhs, Sum&& rhs);
MINIMUM_LINEAR_API Sum operator-(Sum&& lhs, const Sum& rhs);
MINIMUM_LINEAR_API Sum operator-(Sum&& lhs, Sum&& rhs);

namespace {
inline std::ostream& operator<<(std::ostream& out, const Sum& sum) {
	out << sum.value();
	return out;
}
}  // namespace

MINIMUM_LINEAR_API Sum sum(const std::vector<Variable>& xs);
MINIMUM_LINEAR_API Sum sum(const std::vector<BooleanVariable>& xs);

// A sum of several variables (plus a constant).
//
// It is created by the user in these ways:
//
//		auto  l1 = !x;
//		auto  l2 = x || y;
//		auto  l3 = x || !y || !z || w;
//
class MINIMUM_LINEAR_API LogicalExpression {
	friend class Constraint;
	friend class IP;

   public:
	LogicalExpression(const LogicalExpression& expr);
	LogicalExpression& operator=(const LogicalExpression& expr);
	LogicalExpression(LogicalExpression&& expr);
	LogicalExpression& operator=(LogicalExpression&& expr);
	LogicalExpression(const BooleanVariable& variable, bool is_negated = false);

	LogicalExpression& operator|=(const LogicalExpression& lhs);

	const Sum get_sum() const { return sum; }

   private:
	Sum sum;
};

MINIMUM_LINEAR_API LogicalExpression operator||(const LogicalExpression& lhs,
                                                const LogicalExpression& rhs);
MINIMUM_LINEAR_API LogicalExpression operator||(LogicalExpression&& lhs,
                                                const LogicalExpression& rhs);
MINIMUM_LINEAR_API LogicalExpression operator||(const LogicalExpression& lhs,
                                                LogicalExpression&& rhs);
MINIMUM_LINEAR_API LogicalExpression operator||(LogicalExpression&& lhs, LogicalExpression&& rhs);

MINIMUM_LINEAR_API LogicalExpression operator!(const BooleanVariable& variable);

MINIMUM_LINEAR_API LogicalExpression implication(const BooleanVariable& antecedent,
                                                 const LogicalExpression& consequent);
MINIMUM_LINEAR_API LogicalExpression implication(const BooleanVariable& antecedent,
                                                 LogicalExpression&& consequent);
}  // namespace linear
}  // namespace minimum
