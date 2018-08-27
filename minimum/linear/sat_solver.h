#pragma once

#include <memory>
#include <utility>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/linear/ip.h>

namespace minimum {
namespace linear {

class SatSolver {
   public:
	typedef int Var;
	typedef std::pair<Var, bool> Literal;

	virtual ~SatSolver() {}

	// Core methods that must be implemented.
	virtual Var add_variable() = 0;
	virtual void add_clause(const std::vector<Literal>& clause) = 0;
	virtual bool solve(std::vector<int>* solution) = 0;
	virtual bool solve(const std::vector<Literal>& assumptions, std::vector<int>* solution) = 0;

	// Methods below may be overridden if needed.

	// A helper variable can not be used as an assumption later.
	// It may have been eliminated.
	virtual Var add_helper_variable() { return add_variable(); }
	// Simplification can be turned off if the problem is known to be
	// optimally encoded.
	virtual void use_simplification(bool use) {}

	virtual bool could_be_ok(Var var, bool assignment) { return true; }
	virtual bool solve_limited(std::vector<int>* solution, int attempts) { return solve(solution); }
	virtual std::vector<std::vector<std::pair<int, bool>>> get_learnts() { return {}; }

	// Helper methods that are already implemented.
	void add_clause(Literal l1) { add_clause(std::vector<Literal>{l1}); }
	void add_clause(Literal l1, Literal l2) { add_clause(std::vector<Literal>{l1, l2}); }
};

std::unique_ptr<SatSolver> MINIMUM_LINEAR_API minisat_solver(bool silent_mode = false);
std::unique_ptr<SatSolver> MINIMUM_LINEAR_API cleaneling_solver();
std::unique_ptr<SatSolver> MINIMUM_LINEAR_API glucose_solver();
}  // namespace linear
}  // namespace minimum
