#pragma once
#include <cctype>
#include <cmath>
#include <complex>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include <minimum/core/enum.h>
#include <minimum/core/export.h>
#include <minimum/core/variant.h>

namespace minimum::core {

MAKE_ENUM(
    Operation, ADD, SUBTRACT, NEGATE, MULTIPLY, DIVIDE, EXP, LOG, LOG10, POW, SIN, COS, TAN, SQRT);
using Command = std::variant<Operation, double, std::string>;

// Recursive descent parser that can be extended.
//
// Currently parses mathematical expressions into reverse-Polish notation that
// can be easily executed (see evaluate below).
class MINIMUM_CORE_API Parser {
   public:
	Parser(std::string_view data_) : data(data_) {}

	// Parses the expression and checks that no extra tokens remain.
	std::vector<Command> parse();

   protected:
	// Parses expressions of form
	// sum =  product ("+-" product)*
	void parse_sum();

	// Parses expressions of form
	// parenthesized_sum = "(" sum ")"
	// or
	// parenthesized_sum = "(" sum "," sum ")"
	// If more are desired.
	void parse_parenthesized_sum(int num_comma_separated_sums = 1);

	// Parses expressions of the form
	// product = factor ("*" factor)* ("/" factor)?
	void parse_product();

	// Parses expressions of the form
	// term = atom | atom "^" atom
	void parse_factor();

	// Parses expressions of the form
	// atom = "+" atom | "-"" atom | parenthesized_sum | double | symbol
	void parse_atom();

	// Parses a symbol [a-z][a-z0-9]*
	void parse_symbol();

	// Parses a floating-point number.
	void parse_double();

	// Moves past all whitespace.
	void skip_whitespace();

	// Gets the next character or 0 if at the end of the string.
	char get();

	// Like get() but verifies the character.
	void require(char expected, std::string_view error_message);
	void require(char expected);
	void require(std::string_view expected);

	// Returns the next character without advancing.
	char peek() const;

	// Returns true and removes the character if it matches, otherwise
	// does nothing and returns false.
	bool consume(char expected);

	// Returns true and removes the characters if they match, otherwise
	// does nothing and returns false.
	bool consume(std::string_view expected);

	std::string_view data;
	std::vector<Command> result;
};

namespace {
// Converts e.g. "exp" to Operation::EXP.
Operation function_op_from_string(std::string_view s) {
	if (s == "exp") {
		return Operation::EXP;
	} else if (s == "log") {
		return Operation::LOG;
	} else if (s == "log10") {
		return Operation::LOG10;
	} else if (s == "pow") {
		return Operation::POW;
	} else if (s == "sin") {
		return Operation::SIN;
	} else if (s == "cos") {
		return Operation::COS;
	} else if (s == "tan") {
		return Operation::TAN;
	} else if (s == "sqrt") {
		return Operation::SQRT;
	}
	check(false, "Unknown function ", s);
	return {};
}

// How many stack arguments a function opcode (EXP, POW, ...) requires.
constexpr int num_function_arguments(Operation op) {
	switch (op) {
		case Operation::ADD:
		case Operation::SUBTRACT:
		case Operation::MULTIPLY:
		case Operation::DIVIDE:
		case Operation::POW:
			return 2;
		case Operation::NEGATE:
		case Operation::EXP:
		case Operation::LOG:
		case Operation::LOG10:
		case Operation::SIN:
		case Operation::COS:
		case Operation::TAN:
		case Operation::SQRT:
			return 1;
	}
	return -1;
}
}  // namespace

// Stores a list of commands for rapid evaluation.
class Expression {
   public:
	Expression(std::vector<Command> commands_) : commands(std::move(commands_)) {
		// The constructor computes the amount of stack space needed and
		// validates the list of commands. That way, the evaluate function
		// can be implemented without checks and reallocations.
		max_stack_size = 0;
		int stack_size = 0;
		for (auto& command : commands) {
			minimum::core::switch_(
			    command,
			    [&stack_size](double constant) { stack_size += 1; },
			    [&stack_size](Operation op) {
				    int nargs = num_function_arguments(op);
				    minimum_core_assert(nargs >= 0, "Unknown opcode ", to_string(op));
				    minimum_core_assert(stack_size >= nargs, "Stack is too small.");
				    stack_size -= nargs;
				    stack_size += 1;
			    },
			    [&stack_size](const std::string& identifier) { stack_size += 1; });
			max_stack_size = std::max(stack_size, max_stack_size);
		}
		minimum_core_assert(stack_size == 1, "Commands did not evaluate completely.");
	}
	Expression() : commands({Command{0.0}}) {}
	Expression(const Expression&) = default;
	Expression& operator=(const Expression&) = default;

	// Evaluates a list of commands on an initially empty stack.
	// The number type is a template parameter, so can be used for e.g. automatic differentiation.
	template <typename T>
	T evaluate(const std::unordered_map<std::string, T>& identifiers = {}) const {
		using std::cos;
		using std::exp;
		using std::log;
		using std::pow;
		using std::sin;
		using std::sqrt;
		using std::tan;
		const double log10 = std::log(10.0);

		constexpr int static_stack_size = 10;
		T static_stack[static_stack_size + 1];
		T* stack = static_stack;
		std::vector<T> dynamic_stack;
		if (max_stack_size > static_stack_size) {
			dynamic_stack.resize(max_stack_size + 1);
			stack = dynamic_stack.data();
		}
		auto pop = [&stack]() -> T { return *(stack--); };
		auto push = [&stack](T t) -> void { *(++stack) = std::move(t); };

		for (auto& command : commands) {
			minimum::core::switch_(
			    command,
			    [&](double constant) { push(constant); },
			    [&](Operation op) {
				    switch (op) {
					    case Operation::ADD:
					    case Operation::SUBTRACT:
					    case Operation::MULTIPLY:
					    case Operation::DIVIDE:
					    case Operation::POW: {
						    auto right = pop();
						    auto left = pop();
						    if (op == Operation::ADD) {
							    push(left + right);
						    } else if (op == Operation::SUBTRACT) {
							    push(left - right);
						    } else if (op == Operation::MULTIPLY) {
							    push(left * right);
						    } else if (op == Operation::DIVIDE) {
							    push(left / right);
						    } else if (op == Operation::POW) {
							    push(pow(left, right));
						    }
						    break;
					    }
					    case Operation::NEGATE:
						    *stack = -(*stack);
						    break;
					    case Operation::EXP:
						    *stack = exp(*stack);
						    break;
					    case Operation::LOG:
						    *stack = log(*stack);
						    break;
					    case Operation::LOG10:
						    *stack = log(*stack) / log10;
						    break;
					    case Operation::SIN:
						    *stack = sin(*stack);
						    break;
					    case Operation::COS:
						    *stack = cos(*stack);
						    break;
					    case Operation::TAN:
						    *stack = tan(*stack);
						    break;
					    case Operation::SQRT:
						    *stack = sqrt(*stack);
						    break;
					    default:
						    minimum_core_assert(false, "Unknown opcode ", to_string(op));
				    }
			    },
			    [&](const std::string& identifier) {
				    auto itr = identifiers.find(identifier);
				    check(itr != identifiers.end(), "Unknown identifier: ", identifier);
				    push(itr->second);
			    });
		}
		return pop();
	}

   private:
	std::vector<Command> commands;
	int max_stack_size = 0;
};

MINIMUM_CORE_API std::ostream& operator<<(std::ostream& out, const Command& command);
}  // namespace minimum::core
