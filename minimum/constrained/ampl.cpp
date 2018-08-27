#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <streambuf>
#include <string>

#include <minimum/constrained/ampl.h>
#include <minimum/constrained/successive.h>
#include <minimum/core/check.h>
#include <minimum/core/numeric.h>
#include <minimum/core/parser.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/nonlinear/auto_diff_term.h>
#include <minimum/nonlinear/constrained_function.h>
#include <minimum/nonlinear/data/util.h>
#include <minimum/nonlinear/dynamic_auto_diff_term.h>
#include <minimum/nonlinear/solver.h>
#include <minimum/nonlinear/string_functor.h>
using minimum::core::check;
using minimum::core::range;
using minimum::core::to_string;

namespace minimum::constrained {

class LowerBound {
   public:
	LowerBound(int index_, double c_) : index(index_), c{c_} {}

	template <typename R>
	R operator()(const std::vector<int>& dimensions, R* const* const x) const {
		return -x[0][index] + c;
	}

   private:
	int index;
	double c;
};

class UpperBound {
   public:
	UpperBound(int index_, double c_) : index(index_), c{c_} {}

	template <typename R>
	R operator()(const std::vector<int>& dimensions, R* const* const x) const {
		return x[0][index] - c;
	}

   private:
	int index;
	double c;
};

AmplProblem::AmplProblem(std::vector<Variable> variables_,
                         minimum::core::Expression objective_,
                         std::vector<Constraint> constraints_)
    : variables(variables_),
      values(variables_.size()),
      objective(objective_),
      constraints(constraints_) {
	for (auto i : range(variables.size())) {
		values[i] = variables[i].start_value;
		data.emplace_back(&values[i]);
	}
	values_data = values.data();
}

void AmplProblem::reset() {
	for (auto i : range(variables.size())) {
		values[i] = variables[i].start_value;
	}
}

void AmplProblem::create_function(minimum::nonlinear::ConstrainedFunction* function) const {
	using minimum::nonlinear::AutoDiffTerm;
	using minimum::nonlinear::DynamicAutoDiffTerm;
	using minimum::nonlinear::StringFunctor1;

	std::vector<std::string> names;
	const std::vector<int> dimensions = {int(variables.size())};
	for (auto i : range(variables.size())) {
		auto& var = variables[i];

		names.emplace_back(var.name);
		if (var.lower_bound > std::numeric_limits<double>::min()) {
			function->add_constraint_term(
			    to_string(var.name, "-lower"),
			    std::make_shared<DynamicAutoDiffTerm<LowerBound>>(dimensions, i, var.lower_bound),
			    values_data);
		}
		if (var.upper_bound < std::numeric_limits<double>::max()) {
			function->add_constraint_term(
			    to_string(var.name, "-upper"),
			    std::make_shared<DynamicAutoDiffTerm<UpperBound>>(dimensions, i, var.upper_bound),
			    values_data);
		}
	}
	auto objective_term =
	    std::make_shared<DynamicAutoDiffTerm<StringFunctor1>>(dimensions, objective, names);
	function->add_term(objective_term, values_data);

	for (auto& constraint : constraints) {
		auto constraint_term = std::make_shared<DynamicAutoDiffTerm<StringFunctor1>>(
		    dimensions, constraint.expression, names);
		if (constraint.type == minimum::nonlinear::Constraint::Type::EQUAL) {
			function->add_equality_constraint_term(constraint.name, constraint_term, values_data);
		} else if (constraint.type == minimum::nonlinear::Constraint::Type::LESS_THAN_OR_EQUAL) {
			function->add_constraint_term(constraint.name, constraint_term, values_data);
		}
	}
}

void AmplProblem::create_function(minimum::nonlinear::Function* function) const {
	using minimum::nonlinear::AutoDiffTerm;
	using minimum::nonlinear::DynamicAutoDiffTerm;
	using Interval = minimum::nonlinear::Interval<double>;
	using minimum::nonlinear::StringFunctor1;

	minimum_core_assert(constraints.empty());
	function->hessian_is_enabled = false;

	std::vector<std::string> names;
	const std::vector<int> dimensions = {int(variables.size())};
	std::vector<Interval> bounds(variables.size());
	for (auto i : range(variables.size())) {
		auto& var = variables[i];
		names.emplace_back(var.name);

		double lower = -Interval::infinity;
		double upper = Interval::infinity;
		if (var.lower_bound > std::numeric_limits<double>::min()) {
			lower = var.lower_bound;
		}
		if (var.upper_bound < std::numeric_limits<double>::max()) {
			upper = var.upper_bound;
		}
		bounds[i] = {lower, upper};
	}
	function->add_variable(values_data, variables.size());
	function->set_variable_bounds(values_data, bounds);
	auto objective_term =
	    std::make_shared<DynamicAutoDiffTerm<StringFunctor1>>(dimensions, objective, names);
	function->add_term(objective_term, values_data);
}

AmplProblem AmplParser::parse_ampl() {
	try {
		parse_vars();
		parse_objective();
		parse_constraints();
		skip_whitespace();
		check(data.empty(), "Data not completely parsed: ", data);
	} catch (std::exception& error) {
		std::string message = error.what();
		if (message.size() > 100) {
			message = message.substr(0, 100);
		}
		check(false, "Error on line ", current_line(), ": ", message);
	}
	return AmplProblem(std::move(variables), std::move(objective), std::move(constraints));
}

void AmplParser::skip_comments() {
	skip_whitespace();
	while (consume('#')) {
		while (true) {
			char ch = get();
			if (ch == '\n' || ch == '\r') {
				break;
			}
		}
		skip_whitespace();
	}
}

void AmplParser::parse_vars() {
	skip_comments();
	while (consume("var")) {
		AmplProblem::Variable var;
		parse_symbol();
		var.name = std::get<std::string>(result.back());

		skip_whitespace();
		if (consume(":=")) {
			skip_whitespace();
			parse_double();
			var.start_value = std::get<double>(result.back());
			skip_whitespace();
			consume(",");
		}
		skip_whitespace();
		if (peek() != ';') {
			do {
				skip_whitespace();
				if (consume(">=")) {
					skip_whitespace();
					parse_double();
					var.lower_bound = std::get<double>(result.back());
				} else if (consume("<=")) {
					skip_whitespace();
					parse_double();
					var.upper_bound = std::get<double>(result.back());
				} else {
					check(false, "Expected variable bounds.");
				}
				skip_whitespace();
			} while (consume(','));
		}

		require(';', "Expected semicolon ");
		skip_comments();
		variables.emplace_back(std::move(var));
	}
	result.clear();
}
void AmplParser::parse_objective() {
	skip_comments();
	bool maximize;
	if (consume("minimize")) {
		maximize = false;
	} else if (consume("maximize")) {
		maximize = true;
	} else {
		// No objective is allowed. The default expression is 0.
		return;
	}
	skip_comments();
	require("obj");
	skip_comments();
	require(":");
	skip_comments();

	parse_sum();
	require(";");
	auto commands = get_commands();
	if (maximize) {
		commands.emplace_back(minimum::core::Operation::NEGATE);
	}
	objective = minimum::core::Expression(commands);
}
void AmplParser::parse_constraints() {
	skip_comments();
	if (!consume("subject to")) {
		return;
	}
	skip_comments();
	while (isalpha(peek())) {
		AmplProblem::Constraint constraint;
		parse_symbol();
		constraint.name = std::get<std::string>(result.back());
		result.clear();

		skip_whitespace();
		require(":");
		skip_whitespace();
		parse_sum();
		auto commands = get_commands();

		skip_whitespace();
		bool negate = false;
		if (consume("<=")) {
			constraint.type = minimum::nonlinear::Constraint::Type::LESS_THAN_OR_EQUAL;
			negate = true;
		} else if (consume(">=")) {
			constraint.type = minimum::nonlinear::Constraint::Type::LESS_THAN_OR_EQUAL;
		} else if (consume("=")) {
			constraint.type = minimum::nonlinear::Constraint::Type::EQUAL;
		} else {
			check(false, "Unknown constraint type: ", data);
		}
		skip_whitespace();
		parse_sum();
		commands.insert(commands.begin(), result.begin(), result.end());
		result.clear();
		commands.emplace_back(minimum::core::Operation::SUBTRACT);
		if (negate) {
			commands.emplace_back(minimum::core::Operation::NEGATE);
		}
		constraint.expression = minimum::core::Expression(std::move(commands));
		constraints.emplace_back(std::move(constraint));

		skip_whitespace();
		require(";");

		skip_comments();
		consume("subject to");
		skip_comments();
	}
}

int AmplParser::current_line() const {
	return std::count(original_data.begin(),
	                  original_data.begin() + (original_data.size() - data.size()),
	                  '\n')
	       + 1;
}

std::vector<minimum::core::Command> AmplParser::get_commands() {
	std::vector<minimum::core::Command> cmds;
	cmds.swap(result);
	return cmds;
}

}  // namespace minimum::constrained