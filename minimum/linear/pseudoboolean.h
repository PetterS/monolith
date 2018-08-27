#pragma once
#include <memory>

#include <minimum/linear/export.h>
#include <minimum/linear/sum.h>
#include <minimum/linear/variable.h>

namespace minimum {
namespace linear {

// A polynomial of serveral boolean variables.
//
// For example:
//    auto f = x*y*z - 5*y*z + x
//
// Since only boolean variables are allowed, x and x*x are represented
// in the same way internally.
class MINIMUM_LINEAR_API PseudoBoolean {
	friend class IP;

   public:
	PseudoBoolean();
	PseudoBoolean(Variable var);
	PseudoBoolean(const Sum& sum);
	PseudoBoolean(const PseudoBoolean& rhs);
	PseudoBoolean(PseudoBoolean&& rhs);
	~PseudoBoolean();

	PseudoBoolean& operator=(const PseudoBoolean& rhs);
	PseudoBoolean& operator=(PseudoBoolean&& rhs);

	PseudoBoolean& operator+=(const Variable& rhs);
	PseudoBoolean& operator+=(const PseudoBoolean& rhs);
	PseudoBoolean& operator-=(const Variable& rhs);
	PseudoBoolean& operator-=(const PseudoBoolean& rhs);
	void negate();
	PseudoBoolean operator-() const&;
	PseudoBoolean&& operator-() &&;
	PseudoBoolean& operator*=(const PseudoBoolean& rhs);
	PseudoBoolean& operator*=(double rhs);

	double value() const;

	// Returns a string representing the polynomial. Used for
	// testing.
	std::string str() const;

   private:
	const std::map<std::vector<int>, double>& indices() const;

	struct Implementation;
	std::unique_ptr<Implementation> impl;
};

MINIMUM_LINEAR_API PseudoBoolean operator+(Variable lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBoolean operator+(PseudoBoolean lhs, Variable rhs);
// PseudoBoolean with itself.
MINIMUM_LINEAR_API PseudoBoolean operator+(const PseudoBoolean& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator+(PseudoBoolean&& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator+(const PseudoBoolean& lhs, PseudoBoolean&& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator+(PseudoBoolean&& lhs, PseudoBoolean&& rhs);
// PseudoBoolean and sums.
MINIMUM_LINEAR_API PseudoBoolean operator+(Sum, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBoolean operator+(PseudoBoolean lhs, Sum rhs);

MINIMUM_LINEAR_API PseudoBoolean operator-(Variable lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBoolean operator-(PseudoBoolean lhs, Variable rhs);
// PseudoBoolean with itself.
MINIMUM_LINEAR_API PseudoBoolean operator-(const PseudoBoolean& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator-(PseudoBoolean&& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator-(const PseudoBoolean& lhs, PseudoBoolean&& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator-(PseudoBoolean&& lhs, PseudoBoolean&& rhs);
// PseudoBoolean and sums.
MINIMUM_LINEAR_API PseudoBoolean operator-(Sum, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBoolean operator-(PseudoBoolean lhs, Sum rhs);

MINIMUM_LINEAR_API PseudoBoolean operator*(Variable lhs, Variable rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(Variable lhs, PseudoBoolean rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(PseudoBoolean lhs, Variable rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(Variable lhs, Sum rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(Sum lhs, Variable rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(Sum lhs, Sum rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(PseudoBoolean lhs, double rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(double lhs, PseudoBoolean rhs);
// PseudoBoolean with itself.
MINIMUM_LINEAR_API PseudoBoolean operator*(const PseudoBoolean& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(const PseudoBoolean& lhs, PseudoBoolean&& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(PseudoBoolean&& lhs, const PseudoBoolean& rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(PseudoBoolean&& lhs, PseudoBoolean&& rhs);
// PseudoBoolean and sums.
MINIMUM_LINEAR_API PseudoBoolean operator*(PseudoBoolean lhs, Sum rhs);
MINIMUM_LINEAR_API PseudoBoolean operator*(Sum lhs, PseudoBoolean rhs);

MINIMUM_LINEAR_API std::ostream& operator<<(std::ostream& out, const PseudoBoolean& pb);
}  // namespace linear
}  // namespace minimum
