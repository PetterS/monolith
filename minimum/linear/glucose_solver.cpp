#include <iterator>
#include <string>

#include <third-party/glucose-syrup/simp/SimpSolver.h>

#include <minimum/core/check.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/string.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>

using minimum::core::check;

namespace minimum {
namespace linear {

class GlucoseSolver : public SatSolver {
   public:
	GlucoseSolver() { solver.verbosity = -1; }

	virtual Var add_variable() override;
	virtual void add_clause(const std::vector<Literal>& clause) override;
	virtual bool solve(std::vector<int>* solution) override;
	virtual bool solve(const std::vector<Literal>& assumptions,
	                   std::vector<int>* solution) override;

   private:
	// SimpSolver is faster on Nurikabe, but hangs the unit tests with output.
	Glucose::Solver solver;
	std::vector<Glucose::Lit> literals;
};

std::unique_ptr<SatSolver> glucose_solver() { return std::make_unique<GlucoseSolver>(); }

SatSolver::Var GlucoseSolver::add_variable() {
	literals.push_back(Glucose::mkLit(solver.newVar()));
	return literals.size() - 1;
}

void GlucoseSolver::add_clause(const std::vector<Literal>& clause) {
	Glucose::vec<Glucose::Lit> minisat_clause;
	for (auto& literal : clause) {
		if (literal.second) {
			minisat_clause.push(literals.at(literal.first));
		} else {
			minisat_clause.push(~literals.at(literal.first));
		}
	}
	solver.addClause(minisat_clause);
}

bool GlucoseSolver::solve(std::vector<int>* solution) {
	bool result = solver.solve();
	if (result) {
		solution->clear();
		for (auto& lit : literals) {
			auto value = solver.modelValue(lit);
			minimum_core_assert(value == Glucose::l_True || value == Glucose::l_False);
			solution->push_back(value == Glucose::l_True ? 1 : 0);
		}
	}
	return result;
}

bool GlucoseSolver::solve(const std::vector<Literal>& assumptions, std::vector<int>* solution) {
	Glucose::vec<Glucose::Lit> minisat_assumptions;
	for (auto& literal : assumptions) {
		if (literal.second) {
			minisat_assumptions.push(literals.at(literal.first));
		} else {
			minisat_assumptions.push(~literals.at(literal.first));
		}
	}
	bool result = solver.solve(minisat_assumptions);
	if (result) {
		solution->clear();
		for (auto& lit : literals) {
			auto value = solver.modelValue(lit);
			minimum_core_assert(value == Glucose::l_True || value == Glucose::l_False);
			solution->push_back(value == Glucose::l_True ? 1 : 0);
		}
	}
	return result;
}
}  // namespace linear
}  // namespace minimum
