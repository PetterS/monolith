#pragma once

#include <cstddef>
#include <iostream>

#include <minimum/linear/export.h>

namespace minimum {
namespace linear {

class IP;

// Represents a variable in the optimization problem. It is
// only created by the IP class.
class MINIMUM_LINEAR_API Variable {
	friend class IP;
	friend class Sum;
	friend class LogicalExpression;
	friend class BooleanVariable;
	friend class PseudoBoolean;

   public:
	// The value of this variable after the solution has been obtained.
	double value() const;

	// Index of the variable. Normally not needed, but can be useful
	// when a unique ID per variable (per problem) is needed.
	size_t get_index() const { return index; }

   private:
	Variable(std::size_t index_, const IP* creator_) : index(index_), creator(creator_) {}

	std::size_t index;
	const IP* creator;

   public:
	// Everything below ONLY exists so that the variable
	// may be stored in an std::map.
	Variable() : index(0xdeaddead), creator(nullptr) {}
};

// Represents a boolean variable in order to enforce
// compile-time restrictions to boolean expressions.
class MINIMUM_LINEAR_API BooleanVariable : public Variable {
	friend class IP;

   public:
	bool bool_value() const;

   private:
	BooleanVariable(const Variable& variable) : Variable(variable) {}

   public:
	// Everything below ONLY exists so that the variable
	// may be stored in an std::map.
	BooleanVariable() {}
};

// Prints the value of the variable after a solution
// has been obtained.
MINIMUM_LINEAR_API std::ostream& operator<<(std::ostream& out, const Variable& variable);
MINIMUM_LINEAR_API std::ostream& operator<<(std::ostream& out, const BooleanVariable& variable);

class MINIMUM_LINEAR_API DualVariable {
	friend class IP;

   public:
	// The value of this dual variable after the solution has been obtained.
	// Will throw if is_available returns false.
	double value() const;

	// Whether this dual varible corresponds to a constraint. Some constraints
	// are optimized away at the time of creation.
	bool is_valid() const;

	// Whether this dual varible is available in the solution. Some variables
	// may be optimized away or not provided by the solver (e.g. SAT).
	bool is_available() const;

	// Index of the variable. Normally not needed, but can be useful
	// when a unique ID per variable (per problem) is needed.
	size_t get_index() const { return index; }

   private:
	DualVariable() : index(0), creator(nullptr) {}
	DualVariable(std::size_t index_, const IP* creator_) : index(index_), creator(creator_) {}

	std::size_t index;
	const IP* creator;
};

template <typename Var>
class Key {
   public:
	Key(Var var_) : var(var_) {}

	bool operator!=(const Key& rhs) const { return var.get_index() != rhs.var.get_index(); }
	bool operator==(const Key& rhs) const { return var.get_index() == rhs.var.get_index(); }
	bool operator<(const Key& rhs) const { return var.get_index() < rhs.var.get_index(); }

	const Var& get_var() const { return var; }

   private:
	Var var;
};
}  // namespace linear
}  // namespace minimum
