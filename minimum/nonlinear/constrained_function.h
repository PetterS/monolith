// Petter Strandmark.
#ifndef MINIMUM_NONLINEAR_CONSTRAINED_FUNCTION_H
#define MINIMUM_NONLINEAR_CONSTRAINED_FUNCTION_H

#include <minimum/nonlinear/export.h>
#include <minimum/nonlinear/function.h>

namespace minimum {
namespace nonlinear {

class Solver;
struct SolverResults;

class MINIMUM_NONLINEAR_API Constraint : public Function {
   public:
	Constraint() { hessian_is_enabled = false; }

	enum class Type { LESS_THAN_OR_EQUAL, EQUAL };
	Type type = Type::LESS_THAN_OR_EQUAL;
	double lambda = 0;
	std::vector<int> global_indices;

	// Not thread safe.
	// TODO: Change API.
	double evaluate_and_cache() const;

	// Return the value that was stored by evaluate_and_cache.
	double get_cached_value() const;

	// Whether the constraint is feasible with respect to
	// the provided tolerance.
	bool is_feasible(double tolerance) const;

	// Returns the cached value times lambda.
	double infeasibility() const;

	// How much the constraint is violated.
	double violation() const;

	// Update lambda. See Nocedal & Wright on constrained
	// optimization.
	void update_lambda(double mu);

   private:
	mutable double cached_value = 0;
};

//
// Builds an augmented lagrangian of the terms added
// to the objective function and constraints.
//
// See Nocedal & Wright, chapter 17.
//
class MINIMUM_NONLINEAR_API ConstrainedFunction {
   public:
	ConstrainedFunction();

	// Disable copying for now.
	ConstrainedFunction(const ConstrainedFunction&) = delete;
	ConstrainedFunction(ConstrainedFunction&&) = delete;

	~ConstrainedFunction();

	// Function improvement tolerance. The constrained solver stops
	// if |df| / (|f| + tol) < tol
	// between outer iterations.
	double function_improvement_tolerance = 0;

	// When updating the dual variables, the constrained solver
	// stops if the change in lagrangian value divided by the
	// dual variable change is less than this tolerance.
	// (and the feasibility tolerance is satisfied)
	double dual_change_tolerance = 1e-6;

	// Feasibility tolerance. A point is considered feasible if
	// maximum violation is less than this value.
	double feasibility_tolerance = 1e-10;

	// The maximum number of (outer) iterations the solver will
	// perform.
	int max_number_of_iterations = 100;

	// Access to the objective function for e.g. evaluation.
	const Function& objective() const;

	// Function containing all objective and constraint terms. Used
	// by solvers to keep track of indices.
	const Function& all() const;

	// Access to all constraints of the function.
	const std::map<std::string, Constraint>& constraints() const;

	// Whether the provided point is feasible.
	bool is_feasible() const;

	//
	// Adds a term to the objective function.
	//
	void add_term(std::shared_ptr<const Term> term, const std::vector<double*>& arguments);

	template <typename... PointerToDouble>
	void add_term(std::shared_ptr<const Term> term, PointerToDouble... args) {
		add_term(term, {args...});
	}

	template <typename MyTerm, typename... PointerToDouble>
	void add_term(PointerToDouble... args) {
		add_term(std::make_shared<MyTerm>(), {args...});
	}

	//
	// Adds a term as a constraint.
	//
	// This will add the constraint c(x) ≤ 0.
	//
	void add_constraint_term(const std::string& constraint_name,
	                         std::shared_ptr<const Term> term,
	                         const std::vector<double*>& arguments);

	template <typename... PointerToDouble>
	void add_constraint_term(const std::string& constraint_name,
	                         std::shared_ptr<const Term> term,
	                         PointerToDouble... args) {
		add_constraint_term(constraint_name, term, {args...});
	}

	template <typename MyTerm, typename... PointerToDouble>
	void add_constraint_term(const std::string& constraint_name, PointerToDouble... args) {
		add_constraint_term(constraint_name, std::make_shared<MyTerm>(), {args...});
	}

	//
	// Adds a term as a constraint.
	//
	// This will add the constraint c(x) = 0.
	//
	void add_equality_constraint_term(const std::string& constraint_name,
	                                  std::shared_ptr<const Term> term,
	                                  const std::vector<double*>& arguments);

	template <typename... PointerToDouble>
	void add_equality_constraint_term(const std::string& constraint_name,
	                                  std::shared_ptr<const Term> term,
	                                  PointerToDouble... args) {
		add_equality_constraint_term(constraint_name, term, {args...});
	}

	template <typename MyTerm, typename... PointerToDouble>
	void add_equality_constraint_term(const std::string& constraint_name, PointerToDouble... args) {
		add_equality_constraint_term(constraint_name, std::make_shared<MyTerm>(), {args...});
	}

	// Minimized the constrained function. In the future, this method
	// may be extraced to e.g. a ConstrainedSolver class.
	void solve(const Solver& solver, SolverResults* results);

   private:
	void add_constraint_term(const std::string& constraint_name,
	                         Constraint::Type type,
	                         std::shared_ptr<const Term> term,
	                         const std::vector<double*>& arguments);

	class Implementation;
	// unique_pointer would have been nice, but there are issues
	// with sharing these objects across DLL boundaries in VC++.
	Implementation* impl;
};
}  // namespace nonlinear
}  // namespace minimum

#endif
