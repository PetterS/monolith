#pragma once
#include <string>

#include <minimum/constrained/export.h>
#include <minimum/core/parser.h>
#include <minimum/nonlinear/constrained_function.h>

namespace minimum::constrained {

// Holds an AMPL (.mod) problem.
//
// Can create a ConstrainedFunction after parsing.
class MINIMUM_CONSTRAINED_API AmplProblem {
   public:
	struct Variable {
		std::string name;
		double start_value = 0;
		double lower_bound = std::numeric_limits<double>::min();
		double upper_bound = std::numeric_limits<double>::max();
	};

	struct Constraint {
		std::string name;
		minimum::core::Expression expression;
		minimum::nonlinear::Constraint::Type type =
		    minimum::nonlinear::Constraint::Type::LESS_THAN_OR_EQUAL;
	};

	AmplProblem(std::vector<Variable> variables_,
	            minimum::core::Expression objective_,
	            std::vector<Constraint> constraints_);

	const std::vector<Variable>& get_variables() const { return variables; }
	const std::vector<double>& get_values() const { return values; }
	bool has_constraints() const { return !constraints.empty(); }

	// Resets all variables to their original values.
	void reset();

	// Creates a ConstrainedFunction from the data.
	void create_function(minimum::nonlinear::ConstrainedFunction* function) const;

	// Creates a Function from the data if there are no constraints.
	void create_function(minimum::nonlinear::Function* function) const;

   private:
	std::vector<Variable> variables;
	std::vector<double> values;
	std::vector<double*> data;
	double* values_data;
	minimum::core::Expression objective;
	std::vector<Constraint> constraints;
};

// Parses data in AMPL (.mod) format.
class MINIMUM_CONSTRAINED_API AmplParser : private minimum::core::Parser {
   public:
	AmplParser(std::string_view data) : Parser(data), original_data(data) {}
	AmplProblem parse_ampl();

   private:
	void skip_comments();

	void parse_vars();
	void parse_objective();
	void parse_constraints();

	int current_line() const;
	std::vector<minimum::core::Command> get_commands();

	std::string_view original_data;
	std::vector<AmplProblem::Variable> variables;
	minimum::core::Expression objective;
	std::vector<AmplProblem::Constraint> constraints;
};

}  // namespace minimum::constrained