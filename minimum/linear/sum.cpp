
#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/sum.h>
using minimum::core::check;

namespace minimum {
namespace linear {

class Sum::Implementation {
   public:
	Implementation() : constant(0.0), creator(nullptr) {}

	Implementation(const Implementation& rhs) { *this = rhs; }

	Implementation& operator=(const Implementation& rhs) {
		cols = rhs.cols;
		values = rhs.values;
		constant = rhs.constant;
		creator = rhs.creator;
		// std::cerr << "Copy sum.\n";
		// throw std::runtime_error("Copying took place.");
		return *this;
	}

	Implementation(double constant_) : constant(constant_), creator(nullptr) {}

	Implementation(const Variable& variable) : constant(0.0), creator(variable.creator) {
		::minimum::core::check(variable.creator != nullptr,
		                       "Variables used in sums must be created by an IP object.");
		cols.push_back(int(variable.index));
		values.push_back(1.0);
	}

	double constant;
	std::vector<int> cols;
	std::vector<double> values;

	const IP* creator;
};

Sum::Sum() : impl(new Implementation) {}

Sum::Sum(const Sum& sum) : impl(new Implementation(*sum.impl)) {
	// std::clog << "Copy.\n";
}

Sum& Sum::operator=(const Sum& sum) {
	*impl = *sum.impl;
	return *this;
}

Sum::Sum(Sum&& sum) : impl(nullptr) {
	impl = sum.impl;
	sum.impl = nullptr;
	// std::clog << "Move.\n";
}

Sum& Sum::operator=(Sum&& sum) {
	if (impl) {
		delete impl;
	}

	impl = sum.impl;
	sum.impl = nullptr;
	return *this;
}

Sum::Sum(double constant_) : impl(new Implementation(constant_)) {}

Sum::Sum(const Variable& variable) : impl(new Implementation(variable)) {}

Sum::~Sum() {
	if (impl) {
		delete impl;
	}
}

double Sum::value() const {
	if (!impl->creator) {
		// This happens if Sum is constant. No variables
		// have been added.
		minimum_core_assert(impl->cols.empty());
		return impl->constant;
	}

	return impl->creator->get_solution(*this);
}

Sum& Sum::operator+=(const Sum& rhs) {
	match_solvers(rhs);

	impl->constant += rhs.impl->constant;
	for (size_t i = 0; i < rhs.impl->cols.size(); ++i) {
		impl->cols.push_back(rhs.impl->cols[i]);
		impl->values.push_back(rhs.impl->values[i]);
	}
	return *this;
}

Sum& Sum::operator-=(const Sum& rhs) {
	match_solvers(rhs);

	impl->constant -= rhs.impl->constant;
	for (size_t i = 0; i < rhs.impl->cols.size(); ++i) {
		impl->cols.push_back(rhs.impl->cols[i]);
		impl->values.push_back(-rhs.impl->values[i]);
	}
	return *this;
}

Sum& Sum::operator*=(double coeff) {
	if (coeff == 0.0) {
		impl->values.clear();
		impl->cols.clear();
		impl->constant = 0.0;
		return *this;
	}

	for (auto& value : impl->values) {
		value *= coeff;
	}
	impl->constant *= coeff;
	return *this;
}

Sum& Sum::operator/=(double coeff) {
	if (coeff == 0.0) {
		throw std::runtime_error("Sum: Division by zero.");
	}

	*this *= (1.0 / coeff);
	return *this;
}

Sum operator*(double coeff, const Variable& variable) {
	Sum sum(variable);
	sum *= coeff;
	return sum;
}

Sum operator*(const Variable& variable, double coeff) {
	Sum sum(variable);
	sum *= coeff;
	return sum;
}

Sum operator*(double coeff, Sum sum) {
	sum *= coeff;
	return sum;
}

Sum operator*(Sum sum, double coeff) {
	sum *= coeff;
	return sum;
}

Sum operator/(Sum sum, double coeff) {
	sum /= coeff;
	return sum;
}

Sum operator+(Variable lhs, Variable rhs) {
	Sum sum = lhs;
	sum += rhs;
	return sum;
}

Sum operator+(Variable lhs, Sum rhs) {
	rhs += lhs;
	return rhs;
}

Sum operator+(Sum lhs, Variable rhs) {
	lhs += rhs;
	return lhs;
}

Sum operator+(const Sum& lhs, const Sum& rhs) {
	// Copy needed.
	Sum sum(rhs);
	sum += lhs;
	return sum;
}

Sum operator+(const Sum& lhs, Sum&& rhs) {
	Sum sum(std::move(rhs));
	sum += lhs;
	return sum;
}

Sum operator+(Sum&& lhs, const Sum& rhs) {
	Sum sum(std::move(lhs));
	sum += rhs;
	return sum;
}

Sum operator+(Sum&& lhs, Sum&& rhs) {
	Sum sum(std::move(lhs));
	sum += rhs;
	return sum;
}

Sum operator-(Variable lhs, Variable rhs) {
	Sum sum = lhs;
	sum -= rhs;
	return sum;
}

Sum operator-(Variable lhs, Sum rhs) {
	rhs -= lhs;
	rhs.negate();
	return rhs;
}

Sum operator-(Sum lhs, Variable rhs) {
	lhs -= rhs;
	return lhs;
}

Sum operator-(Sum sum) {
	sum.negate();
	return sum;
}

Sum operator-(const Variable& variable) {
	Sum sum(variable);
	sum.negate();
	return sum;
}

void Sum::negate() {
	impl->constant = -impl->constant;
	for (auto& value : impl->values) {
		value = -value;
	}
}

Sum operator-(const Sum& lhs, const Sum& rhs) {
	// Copy needed.
	Sum sum(lhs);
	sum -= rhs;
	return sum;
}

Sum operator-(const Sum& lhs, Sum&& rhs) {
	Sum sum(std::move(rhs));
	sum.negate();
	sum += lhs;
	return sum;
}

Sum operator-(Sum&& lhs, const Sum& rhs) {
	Sum sum(std::move(lhs));
	sum -= rhs;
	return sum;
}

Sum operator-(Sum&& lhs, Sum&& rhs) {
	Sum sum(std::move(rhs));
	sum.negate();
	sum += lhs;
	return sum;
}

std::size_t Sum::size() const { return impl->cols.size(); }

const int* Sum::indices() const {
	minimum_core_assert(!impl->cols.empty());
	return impl->cols.data();
}
const double* Sum::values() const {
	minimum_core_assert(!impl->values.empty());
	return impl->values.data();
}
double Sum::constant() const { return impl->constant; }
const IP* Sum::creator() const { return impl->creator; }

void Sum::match_solvers(const Sum& sum) {
	check(impl->creator == nullptr || sum.impl->creator == nullptr
	          || impl->creator == sum.impl->creator,
	      "Variables from different solver can not be mixed.");
	if (impl->creator == nullptr) {
		impl->creator = sum.impl->creator;
	}
}

Sum sum(const std::vector<BooleanVariable>& xs) {
	Sum result = 0;
	for (auto& x : xs) {
		result += x;
	}
	return result;
}

Sum sum(const std::vector<Variable>& xs) {
	Sum result = 0;
	for (auto& x : xs) {
		result += x;
	}
	return result;
}

LogicalExpression operator!(const BooleanVariable& variable) {
	// TODO: check that variable is integer through creator.
	return LogicalExpression(variable, true);
}

LogicalExpression::LogicalExpression(const BooleanVariable& variable, bool is_negated) {
	if (is_negated) {
		sum += 1;
		sum -= variable;
	} else {
		sum += variable;
	}
}

LogicalExpression::LogicalExpression(const LogicalExpression& expr) { *this = expr; }

LogicalExpression& LogicalExpression::operator=(const LogicalExpression& expr) {
	sum = expr.sum;
	return *this;
}

LogicalExpression::LogicalExpression(LogicalExpression&& expr) { *this = std::move(expr); }

LogicalExpression& LogicalExpression::operator=(LogicalExpression&& expr) {
	sum = std::move(expr.sum);
	return *this;
}

LogicalExpression& LogicalExpression::operator|=(const LogicalExpression& lhs) {
	sum.match_solvers(lhs.sum);

	sum += lhs.sum;
	return *this;
}

LogicalExpression operator||(const LogicalExpression& lhs, const LogicalExpression& rhs) {
	// Copying necessary.
	LogicalExpression result(lhs);
	result |= rhs;
	return result;
}

LogicalExpression operator||(LogicalExpression&& lhs, const LogicalExpression& rhs) {
	LogicalExpression result(std::move(lhs));
	result |= rhs;
	return result;
}

LogicalExpression operator||(const LogicalExpression& lhs, LogicalExpression&& rhs) {
	LogicalExpression result(std::move(rhs));
	result |= lhs;
	return result;
}

LogicalExpression operator||(LogicalExpression&& lhs, LogicalExpression&& rhs) {
	LogicalExpression result(std::move(lhs));
	result |= rhs;
	return result;
}

LogicalExpression implication(const BooleanVariable& antecedent,
                              const LogicalExpression& consequent) {
	return !antecedent || consequent;
}

LogicalExpression implication(const BooleanVariable& antecedent, LogicalExpression&& consequent) {
	return !antecedent || std::move(consequent);
}
}  // namespace linear
}  // namespace minimum
