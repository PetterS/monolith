#include <cctype>
#include <cmath>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/parser.h>
#include <minimum/core/variant.h>
using minimum::core::check;
using minimum::core::switch_;

namespace minimum::core {

std::vector<Command> Parser::parse() {
	parse_sum();
	skip_whitespace();
	check(data.empty(), "Syntax error: ", data);

	std::vector<Command> result_to_return;
	result_to_return.swap(result);
	return result_to_return;
}

void Parser::parse_sum() {
	parse_product();
	skip_whitespace();

	while (peek() == '+' || peek() == '-') {
		auto op = get();
		parse_product();
		skip_whitespace();
		if (op == '+') {
			result.emplace_back(Operation::ADD);
		} else if (op == '-') {
			result.emplace_back(Operation::SUBTRACT);
		}
	}
}

void Parser::parse_parenthesized_sum(int num_comma_separated_sums) {
	skip_whitespace();
	minimum_core_assert(get() == '(');

	for (int i = 0; i < num_comma_separated_sums; ++i) {
		parse_sum();
		skip_whitespace();
		if (i < num_comma_separated_sums - 1) {
			require(',', "Expected comma.");
		}
	}
	require(')', "Missing closing bracket.");
}

void Parser::parse_product() {
	parse_factor();
	skip_whitespace();

	while (consume('*')) {
		parse_factor();
		skip_whitespace();
		result.emplace_back(Operation::MULTIPLY);
	}
	if (consume('/')) {
		parse_factor();
		skip_whitespace();
		result.emplace_back(Operation::DIVIDE);
	}
}

void Parser::parse_factor() {
	skip_whitespace();
	parse_atom();
	skip_whitespace();
	if (consume('^') || consume("**")) {
		parse_atom();
		result.emplace_back(Operation::POW);
	}
}

void Parser::parse_atom() {
	skip_whitespace();
	if (peek() == '+' || peek() == '-') {
		auto op = get();
		parse_atom();
		if (op == '-') {
			result.emplace_back(Operation::NEGATE);
		}
	} else if (peek() == '(') {
		parse_parenthesized_sum();
	} else if (std::isdigit(peek()) || (peek() == '.')) {
		parse_double();
	} else if (std::isalpha(peek())) {
		parse_symbol();
	} else {
		check(false, "Expected a number, symbol or parenthesized expression: ", data);
	}
}

void Parser::parse_symbol() {
	std::string symbol;
	skip_whitespace();
	minimum_core_assert(isalpha(peek()), "Expected letter at beginning of symbol: ", data);
	while (isalnum(peek())) {
		symbol.push_back(get());
	}

	if (peek() == '(') {
		// This is a function call. Convert it to an opcode.
		auto op = function_op_from_string(symbol);
		int num_args = num_function_arguments(op);
		parse_parenthesized_sum(num_args);
		result.emplace_back(op);
	} else {
		// No function call; just push the variable.
		result.emplace_back(std::move(symbol));
	}
}

void Parser::parse_double() {
	skip_whitespace();
	std::string s;

	while (peek() == '-' || peek() == '+') {
		s += get();
	}

	// Before the decimal separator.
	while (std::isdigit(peek())) {
		s += get();
	}
	// After the decimal separator.
	if (peek() == '.') {
		s.push_back(get());
		while (std::isdigit(peek())) {
			s.push_back(get());
		}
	}
	// Possible exponent.
	if (peek() == 'e' || peek() == 'E') {
		s.push_back(get());
		if (peek() == '+' || peek() == '-') {
			s.push_back(get());
		}
		check(std::isdigit(peek()), "Expected exponent in floating point constant.");
		while (std::isdigit(peek())) {
			s.push_back(get());
		}
	}
	double number;
	try {
		number = std::stod(s);
	} catch (std::exception&) {
		check(false, s, " is not a number.");
	}
	result.emplace_back(number);
}

void Parser::skip_whitespace() {
	data.remove_prefix(std::min(data.find_first_not_of(" \t\n\r"), data.size()));
}

char Parser::get() {
	if (data.empty()) {
		return 0;
	} else {
		char result = data.front();
		data = data.substr(1);
		return result;
	}
}

void Parser::require(char expected, std::string_view error_message) {
	char result = get();
	check(result == expected, error_message, result);
}

void Parser::require(char expected) {
	char result = get();
	check(result == expected, "Expected ", expected, "; got ", result);
}

void Parser::require(std::string_view expected) {
	check(consume(expected), "Expected ", expected, " at ", data);
}

char Parser::peek() const {
	if (data.empty()) {
		return 0;
	} else {
		return data[0];
	}
}

bool Parser::consume(char expected) {
	if (peek() == expected) {
		get();
		return true;
	}
	return false;
}

bool Parser::consume(std::string_view expected) {
	if (data.size() >= expected.size() && data.compare(0, expected.size(), expected) == 0) {
		data = data.substr(expected.size());
		return true;
	}
	return false;
}

std::ostream& operator<<(std::ostream& out, const Command& command) {
	switch_(command,
	        [&out](double constant) { out << constant; },
	        [&out](Operation op) { out << to_string(op); },
	        [&out](const std::string& identifier) { out << identifier; });

	return out;
}

}  // namespace minimum::core