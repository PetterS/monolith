
#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/pseudoboolean.h>
using minimum::core::check;

namespace minimum {
namespace linear {

struct PseudoBoolean::Implementation {
	Implementation() : creator(nullptr) {}

	void remove_zeros() const;
	mutable std::map<std::vector<int>, double> indices;
	const IP* creator;

	void match_solvers(const Variable& var);
	void match_solvers(const PseudoBoolean& pb);
	void check_type(const Variable& variable) const;
};

void PseudoBoolean::Implementation::match_solvers(const Variable& rhs) {
	check(creator == nullptr || rhs.creator == nullptr || creator == rhs.creator,
	      "Variables from different solvers can not be mixed.");
	if (creator == nullptr) {
		creator = rhs.creator;
	}
}

void PseudoBoolean::Implementation::check_type(const Variable& variable) const {
	minimum_core_assert(variable.creator);
	minimum_core_assert(variable.index < variable.creator->get().variable_size());
	auto& var = variable.creator->get().variable(variable.index);
	check(var.type() == proto::Variable_Type_INTEGER
	          && (var.bound().lower() == 0 || var.bound().lower() == 1)
	          && (var.bound().upper() == 0 || var.bound().upper() == 1),
	      "Variables in a PseudoBoolean expression need to be in {0, 1}.");
}

void PseudoBoolean::Implementation::match_solvers(const PseudoBoolean& rhs) {
	check(creator == nullptr || rhs.impl->creator == nullptr || creator == rhs.impl->creator,
	      "Variables from different solvers can not be mixed.");
	if (creator == nullptr) {
		creator = rhs.impl->creator;
	}
}

PseudoBoolean::PseudoBoolean() : impl(new Implementation()) {}

PseudoBoolean::PseudoBoolean(Variable var) : impl(new Implementation()) {
	impl->check_type(var);
	impl->creator = var.creator;
	impl->indices[{(int)var.index}] = 1;
}

PseudoBoolean::PseudoBoolean(const Sum& sum) : impl(new Implementation()) {
	impl->creator = sum.creator();
	if (sum.constant() != 0) {
		impl->indices[{}] += sum.constant();
	}
	for (size_t i = 0; i < sum.size(); ++i) {
		auto index = sum.indices()[i];
		impl->check_type(Variable(index, sum.creator()));
		impl->indices[{index}] += sum.values()[i];
	}
}

PseudoBoolean::PseudoBoolean(const PseudoBoolean& rhs) : impl(new Implementation()) {
	impl->indices = rhs.impl->indices;
	impl->creator = rhs.impl->creator;
	// std::clog << "Copy.\n";
}

PseudoBoolean::PseudoBoolean(PseudoBoolean&& rhs) : impl(std::move(rhs.impl)) {
	// std::clog << "Move.\n";
}

PseudoBoolean::~PseudoBoolean() {}

PseudoBoolean& PseudoBoolean::operator=(const PseudoBoolean& rhs) {
	impl->indices = rhs.impl->indices;
	impl->creator = rhs.impl->creator;
	// std::clog << "Copy.\n";
	return *this;
}

PseudoBoolean& PseudoBoolean::operator=(PseudoBoolean&& rhs) {
	impl = std::move(rhs.impl);
	return *this;
}

PseudoBoolean& PseudoBoolean::operator+=(const Variable& rhs) {
	impl->check_type(rhs);
	impl->match_solvers(rhs);
	impl->indices[{(int)rhs.index}] += 1;
	return *this;
}

PseudoBoolean& PseudoBoolean::operator+=(const PseudoBoolean& rhs) {
	impl->match_solvers(rhs);
	for (auto& entry : rhs.impl->indices) {
		impl->indices[entry.first] += entry.second;
	}
	return *this;
}

PseudoBoolean operator+(Variable lhs, PseudoBoolean rhs) {
	rhs += lhs;
	return rhs;
}

PseudoBoolean operator+(PseudoBoolean lhs, Variable rhs) {
	lhs += rhs;
	return lhs;
}

PseudoBoolean operator+(const PseudoBoolean& lhs, const PseudoBoolean& rhs) {
	PseudoBoolean result = lhs;
	result += rhs;
	// std::clog << "operator + (const PseudoBoolean & lhs, const PseudoBoolean & rhs)\n";
	return result;
}

PseudoBoolean operator+(PseudoBoolean&& lhs, const PseudoBoolean& rhs) {
	lhs += rhs;
	return std::move(lhs);
}

PseudoBoolean operator+(const PseudoBoolean& lhs, PseudoBoolean&& rhs) {
	rhs += lhs;
	return std::move(rhs);
}

PseudoBoolean operator+(PseudoBoolean&& lhs, PseudoBoolean&& rhs) {
	lhs += rhs;
	return std::move(lhs);
}

PseudoBoolean operator+(Sum lhs, PseudoBoolean rhs) {
	rhs += std::move(lhs);
	return rhs;
}

PseudoBoolean operator+(PseudoBoolean lhs, Sum rhs) {
	lhs += std::move(rhs);
	return lhs;
}

PseudoBoolean& PseudoBoolean::operator-=(const Variable& rhs) {
	impl->check_type(rhs);
	impl->match_solvers(rhs);
	impl->indices[{(int)rhs.index}] -= 1;
	return *this;
}

PseudoBoolean& PseudoBoolean::operator-=(const PseudoBoolean& rhs) {
	impl->match_solvers(rhs);
	for (auto& entry : rhs.impl->indices) {
		impl->indices[entry.first] -= entry.second;
	}
	return *this;
}

PseudoBoolean operator-(Variable lhs, PseudoBoolean rhs) {
	rhs -= lhs;
	rhs.negate();
	return rhs;
}

PseudoBoolean operator-(PseudoBoolean lhs, Variable rhs) {
	lhs -= rhs;
	return lhs;
}

PseudoBoolean operator-(const PseudoBoolean& lhs, const PseudoBoolean& rhs) {
	PseudoBoolean result = lhs;
	result -= rhs;
	return result;
}

PseudoBoolean operator-(PseudoBoolean&& lhs, const PseudoBoolean& rhs) {
	lhs -= rhs;
	return std::move(lhs);
}

PseudoBoolean operator-(const PseudoBoolean& lhs, PseudoBoolean&& rhs) {
	rhs -= lhs;
	rhs.negate();
	return std::move(rhs);
}

PseudoBoolean operator-(PseudoBoolean&& lhs, PseudoBoolean&& rhs) {
	lhs -= rhs;
	return std::move(lhs);
}

PseudoBoolean operator-(Sum lhs, PseudoBoolean rhs) {
	rhs -= std::move(lhs);
	rhs.negate();
	return rhs;
}

PseudoBoolean operator-(PseudoBoolean lhs, Sum rhs) {
	lhs -= std::move(rhs);
	return lhs;
}

void PseudoBoolean::negate() {
	for (auto& entry : impl->indices) {
		entry.second = -entry.second;
	}
}

PseudoBoolean PseudoBoolean::operator-() const& {
	PseudoBoolean result = *this;
	result.negate();
	return result;
}

PseudoBoolean&& PseudoBoolean::operator-() && {
	negate();
	return std::move(*this);
}

PseudoBoolean& PseudoBoolean::operator*=(const PseudoBoolean& rhs) {
	impl->match_solvers(rhs);

	std::map<std::vector<int>, double> indices;
	for (auto& entry : impl->indices) {
		for (auto& rhs_entry : rhs.impl->indices) {
			double coef = entry.second * rhs_entry.second;
			std::vector<int> monomial;
			for (auto& i : entry.first) {
				monomial.push_back(i);
			}
			for (auto& i : rhs_entry.first) {
				monomial.push_back(i);
			}

			std::sort(begin(monomial), end(monomial));
			// Since x*x = x for pseudo-boolean expressions, remove the duplicates.
			monomial.erase(std::unique(begin(monomial), end(monomial)), end(monomial));
			indices[std::move(monomial)] += coef;
		}
	}
	impl->indices = std::move(indices);
	return *this;
}

PseudoBoolean& PseudoBoolean::operator*=(double rhs) {
	for (auto& entry : impl->indices) {
		entry.second *= rhs;
	}
	return *this;
}

double PseudoBoolean::value() const {
	double value = 0;
	for (auto& entry : impl->indices) {
		double monomial_value = 1;
		for (auto index : entry.first) {
			minimum_core_assert(impl->creator != nullptr, "Invalid PseudoBoolean creator.");
			monomial_value *= impl->creator->get_solution(Variable(index, impl->creator));
		}
		value += entry.second * monomial_value;
	}
	return value;
}

PseudoBoolean operator*(Variable lhs, Variable rhs) {
	// TODO: This can be made more efficient.
	PseudoBoolean result = lhs;
	result *= rhs;
	return result;
}

PseudoBoolean operator*(Variable lhs, PseudoBoolean rhs) {
	rhs *= lhs;
	return rhs;
}

PseudoBoolean operator*(PseudoBoolean lhs, Variable rhs) {
	lhs *= rhs;
	return lhs;
}

PseudoBoolean operator*(Variable lhs, Sum rhs) {
	PseudoBoolean result = std::move(rhs);
	result *= lhs;
	return result;
}

PseudoBoolean operator*(Sum lhs, Variable rhs) {
	PseudoBoolean result = std::move(lhs);
	result *= rhs;
	return result;
}

PseudoBoolean operator*(const PseudoBoolean& lhs, const PseudoBoolean& rhs) {
	PseudoBoolean result = lhs;
	result *= rhs;
	// std::clog << "operator * (const PseudoBoolean & lhs, const PseudoBoolean & rhs)\n";
	return result;
}

PseudoBoolean operator*(Sum lhs, Sum rhs) {
	PseudoBoolean result = std::move(lhs);
	result *= std::move(rhs);
	// std::clog << "operator * (const Sum & lhs, const Sum & rhs)\n";
	return result;
}

PseudoBoolean operator*(const PseudoBoolean& lhs, PseudoBoolean&& rhs) {
	rhs *= lhs;
	return std::move(rhs);
}

PseudoBoolean operator*(PseudoBoolean&& lhs, const PseudoBoolean& rhs) {
	lhs *= rhs;
	return std::move(lhs);
}

PseudoBoolean operator*(PseudoBoolean lhs, double rhs) {
	lhs *= rhs;
	return lhs;
}

PseudoBoolean operator*(double lhs, PseudoBoolean rhs) {
	rhs *= lhs;
	return rhs;
}

PseudoBoolean operator*(PseudoBoolean&& lhs, PseudoBoolean&& rhs) {
	lhs *= rhs;
	return std::move(lhs);
}

PseudoBoolean operator*(Sum lhs, PseudoBoolean rhs) {
	rhs *= std::move(lhs);
	return rhs;
}

PseudoBoolean operator*(PseudoBoolean lhs, Sum rhs) {
	lhs *= std::move(rhs);
	return lhs;
}

std::string PseudoBoolean::str() const {
	std::stringstream sout;
	impl->remove_zeros();
	for (auto& entry : impl->indices) {
		if (entry.first.empty()) {
			sout << minimum::core::to_string(entry.second) << " ";
		} else {
			sout << minimum::core::to_string(entry.second) + "*x"
			            + minimum::core::join("*x", entry.first)
			     << " ";
		}
	}
	return sout.str();
}

void PseudoBoolean::Implementation::remove_zeros() const {
	auto itr = indices.begin();
	while (itr != indices.end()) {
		if (itr->second == 0) {
			itr = indices.erase(itr);
		} else {
			itr++;
		}
	}
}

std::ostream& operator<<(std::ostream& out, const PseudoBoolean& pb) {
	check(false, "No printing of PseudoBooleans yet.");
	return out;
}

const std::map<std::vector<int>, double>& PseudoBoolean::indices() const { return impl->indices; }
}  // namespace linear
}  // namespace minimum
