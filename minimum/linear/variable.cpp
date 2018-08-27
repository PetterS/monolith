
#include <minimum/core/check.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/variable.h>

namespace minimum {
namespace linear {

double Variable::value() const { return creator->get_solution(*this); }

bool BooleanVariable::bool_value() const { return creator->get_solution(*this); }

double DualVariable::value() const {
	minimum::core::check(creator != nullptr, "Dual variable not available.");
	return creator->get_solution(*this);
}

bool DualVariable::is_valid() const { return creator != nullptr; }

bool DualVariable::is_available() const {
	if (creator == nullptr) {
		return false;
	}
	auto sol = creator->get_solution(*this);
	return sol == sol;
}

std::ostream& operator<<(std::ostream& out, const Variable& variable) {
	out << variable.value();
	return out;
}

std::ostream& operator<<(std::ostream& out, const BooleanVariable& variable) {
	out << variable.bool_value();
	return out;
}
}  // namespace linear
}  // namespace minimum
