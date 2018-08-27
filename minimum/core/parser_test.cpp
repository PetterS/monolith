#include <cctype>
#include <cmath>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include <catch.hpp>

#include <minimum/core/check.h>

#include <minimum/core/parser.h>
#include <minimum/core/string.h>
#include <minimum/core/variant.h>
using Catch::Contains;
using Catch::Detail::Approx;
using minimum::core::Expression;
using minimum::core::Parser;

namespace {
double evaluate_for_test(std::string_view expr,
                         const std::unordered_map<std::string, double>& identifiers = {}) {
	Expression expression = Parser(expr).parse();
	return expression.evaluate<double>(identifiers);
}
}  // namespace

TEST_CASE("parse_double") {
	CHECK(evaluate_for_test("1.4") == 1.4);
	CHECK(evaluate_for_test("1000") == 1000);
	CHECK_THROWS_WITH(evaluate_for_test("15 + 1.5e + 8"), Contains("Expected exponent"));
}

TEST_CASE("to_string") {
	CHECK(to_string(Parser("1 + 2 * 4").parse()) == "[1, 2, 4, MULTIPLY, ADD]");
}

TEST_CASE("parse_evaluate") {
	CHECK(evaluate_for_test("1 + 2") == 3.0);
	CHECK(evaluate_for_test("1 + 2 * 4") == 9.0);
	CHECK(evaluate_for_test("(1 + (3 * (4 + 5))) + 1000") == 1028);
	CHECK(evaluate_for_test("1 - 5") == -4);
	CHECK(evaluate_for_test("2^3 + (1+1) ^ (0+1+1)") == Approx(8 + 4));
	CHECK(evaluate_for_test("2**3 + (1+1) ** (0+1+1)") == Approx(8 + 4));
	CHECK(evaluate_for_test("1/2") == Approx(0.5));
	CHECK(evaluate_for_test("-1 + (+1)") == 0);
	CHECK(evaluate_for_test("-2^3") == Approx(-8));
	CHECK(evaluate_for_test("2^(-3)") == Approx(0.125));
	CHECK(evaluate_for_test("2^-3") == Approx(0.125));
	CHECK(evaluate_for_test("---3") == -3);
	CHECK(evaluate_for_test("2 + 3") == 5);
	CHECK(evaluate_for_test("2 * 3") == 6);
	CHECK(evaluate_for_test("89") == 89);
	CHECK(evaluate_for_test("   12        -  8   ") == 4);
	CHECK(evaluate_for_test("142        -9   ") == 133);
	CHECK(evaluate_for_test("72+  15") == 87);
	CHECK(evaluate_for_test(" 12*  4") == 48);
	CHECK(evaluate_for_test(" 50/10") == 5);
	CHECK(evaluate_for_test("2.5") == 2.5);
	CHECK(evaluate_for_test("4*2.5 + 8.5+1.5 / 3.0") == 19);
	CHECK(evaluate_for_test("67+2") == 69);
	CHECK(evaluate_for_test(" 2-7") == -5);
	CHECK(evaluate_for_test("5*7 ") == 35);
	CHECK(evaluate_for_test("8/4") == 2);
	CHECK(evaluate_for_test("2 -4") == -2);
	CHECK(evaluate_for_test("2 -4 +6") == 4);
	CHECK(evaluate_for_test("2 -4 +6 -1 -1- 0 +8") == 10);
	CHECK(evaluate_for_test("1 -1   + 2   - 2   +  4 - 4 +    6") == 6);
	CHECK(evaluate_for_test("2 -4 +6 -1 -1- 0 +8") == 10);
	CHECK(evaluate_for_test(" 2*3 - 4*5 + 6/3 ") == -12);
	CHECK(evaluate_for_test("2*2*2*2") == 16);
	CHECK(evaluate_for_test("2*2*2*2/3 + 1") == (16 / 3.0 + 1));
	CHECK(evaluate_for_test("10/4") == 2.5);
	CHECK(evaluate_for_test("(2)") == 2);
	CHECK(evaluate_for_test("(5 + 2*3 - 1 + 7 * 8)") == 66);
	CHECK(evaluate_for_test("(67 + 2 * 3 - 67 + 2/1 - 7)") == 1);
	CHECK(evaluate_for_test("(2) + (17*2-30) * (5)+2 - (8/2)*4") == 8);
	CHECK(evaluate_for_test("(((((5)))))") == 5);
	CHECK(evaluate_for_test("(( ((2)) + 4))*((5))") == 30);
	CHECK(evaluate_for_test("3 + 8/5 -1 -2*5") == Approx(-6.4));
	CHECK(evaluate_for_test("5.0005 + 0.0095") == Approx(5.01));
}

TEST_CASE("large_stack") {
	constexpr int stack_size = 100;
	std::string s;
	for (int i = 0; i < stack_size; ++i) {
		s += "0+(";
	}
	s += "0+0";
	for (int i = 0; i < stack_size; ++i) {
		s += ")";
	}
	CHECK(evaluate_for_test(s) == 0);
}

TEST_CASE("function_call") {
	CHECK(evaluate_for_test("1 + exp(1 + (1 + 1)) + 2") == (1 + std::exp(3.0) + 2));
	CHECK(evaluate_for_test("1 + pow(3, 2)") == Approx(10.0));
	CHECK(evaluate_for_test("exp(0)") == 1);
	CHECK(evaluate_for_test("exp(0)  ") == 1);
	CHECK(evaluate_for_test("exp(0  )") == 1);
	CHECK(evaluate_for_test("exp(  0)") == 1);
	CHECK(evaluate_for_test("sin(0)") == Approx(0));
	CHECK(evaluate_for_test("cos(0)") == Approx(1));
	CHECK(evaluate_for_test("tan(1)") == std::tan(1.0));
}

TEST_CASE("parse_errors") {
	CHECK_THROWS_WITH(Parser("1 + (4 + (5 - 1)").parse(), Contains("Missing closing bracket"));
	CHECK_THROWS_WITH(Parser("1 + ===").parse(),
	                  Contains("Expected a number, symbol or parenthesized expression"));
	CHECK_THROWS_WITH(Parser("1 + 1 petter").parse(), Contains("Syntax error: petter"));
	CHECK_THROWS_WITH(Parser("1 + petter(3)").parse(), Contains("Unknown function petter"));
	CHECK_THROWS_WITH(Parser("1 + pow(3)").parse(), Contains("Expected comma."));

	// Function calls can not have spaces right after the function name.
	CHECK_THROWS_WITH(Parser("exp  (0)").parse(), Contains("Syntax error: (0)"));

	// Disallow ambiguous division.
	CHECK_THROWS_WITH(Parser("5/2*4").parse(), "Syntax error: *4");

	// Bad errors:
	CHECK_THROWS_WITH(Parser("1 + exp(3 = 2)").parse(), Contains("Missing closing bracket"));
	CHECK_THROWS_WITH(Parser("1 + exp(3, 2)").parse(), Contains("Missing closing bracket"));
}

TEST_CASE("evaluate_identifiers") {
	CHECK(evaluate_for_test("x + y", {{"x", 2}, {"y", 1}}) == 3.0);
	CHECK(evaluate_for_test("-x + (+x)", {{"x", 42}}) == Approx(0));
	CHECK_THROWS_WITH(evaluate_for_test("1.0 + a", {{"x", 15}}), Contains("Unknown identifier"));
}

TEST_CASE("complex") {
	// Other numeric types than double should work for the same expression.
	std::unordered_map<std::string, std::complex<double>> vars;
	vars["one"] = {1, 0};
	vars["i"] = {0, 1};
	Expression expr = Parser("one + 2*i").parse();

	CHECK(expr.evaluate<std::complex<double>>(vars) == std::complex<double>(1, 2));

	std::unordered_map<std::string, double> realvars;
	realvars["one"] = 1;
	realvars["i"] = 1;
	CHECK(expr.evaluate<double>(realvars) == 3);
}
