#include <iterator>
#include <string>

#include <third-party/minisat/minisat/core/Solver.h>
#include <third-party/minisat/minisat/simp/SimpSolver.h>

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

class MinisatInternalSatSolver : public SatSolver {
   public:
	MinisatInternalSatSolver(bool silent_mode);

	virtual Var add_variable() override;
	virtual Var add_helper_variable() override;
	virtual void add_clause(const std::vector<Literal>& clause) override;
	virtual bool solve(std::vector<int>* solution) override;
	virtual bool solve(const std::vector<Literal>& assumptions,
	                   std::vector<int>* solution) override;
	virtual void use_simplification(bool use) override { do_simplification = use; }

	// Used to save CNF below.
	Minisat::SimpSolver& get_solver() { return solver; };

	virtual bool could_be_ok(SatSolver::Var var, bool assignment) override;
	virtual std::vector<std::vector<std::pair<int, bool>>> get_learnts() override;

   private:
	Minisat::SimpSolver solver;
	std::vector<Minisat::Lit> literals;
	std::vector<std::vector<std::pair<int, bool>>> learnts;

	bool do_simplification = true;
};

std::unique_ptr<SatSolver> minisat_solver(bool silent_mode) {
	return std::make_unique<MinisatInternalSatSolver>(silent_mode);
}

MinisatInternalSatSolver::MinisatInternalSatSolver(bool silent_mode) {
	solver.verbosity = silent_mode ? 0 : 1;
}

SatSolver::Var MinisatInternalSatSolver::add_variable() {
	auto var = solver.newVar();
	literals.push_back(Minisat::mkLit(var));
	solver.freezeVar(var);
	return literals.size() - 1;
}

SatSolver::Var MinisatInternalSatSolver::add_helper_variable() {
	literals.push_back(Minisat::mkLit(solver.newVar()));
	return literals.size() - 1;
}

bool MinisatInternalSatSolver::could_be_ok(SatSolver::Var var, bool assignment) {
	auto lit = literals.at(var);
	if (!assignment) {
		lit = ~lit;
	}
	Minisat::vec<Minisat::Lit> assumption;
	assumption.push(lit);
	solver.setConfBudget(1);
	solver.verbosity = 0;
	auto result = solver.solveLimited(assumption);
	solver.budgetOff();
	return result != Minisat::l_False;
}

void MinisatInternalSatSolver::add_clause(const std::vector<Literal>& clause) {
	Minisat::vec<Minisat::Lit> minisat_clause;
	for (auto& literal : clause) {
		if (literal.second) {
			minisat_clause.push(literals.at(literal.first));
		} else {
			minisat_clause.push(~literals.at(literal.first));
		}
	}
	solver.addClause(minisat_clause);
}

bool MinisatInternalSatSolver::solve(std::vector<int>* solution) {
	for (int i = 0; i < solution->size(); ++i) {
		solver.setPolarity(Minisat::var(literals.at(i)),
		                   (*solution)[i] > 0 ? Minisat::l_False : Minisat::l_True);
	}

	bool result = solver.solve(do_simplification);
	// solver.learntsize_factor = 1;
	// solver.setConfBudget(100000);
	// bool result = solver.solveLimited(Minisat::vec<Minisat::Lit>()) != Minisat::l_Undef;
	if (result) {
		solution->clear();
		for (auto& lit : literals) {
			auto value = solver.modelValue(lit);
			minimum_core_assert(value == Minisat::l_True || value == Minisat::l_False);
			solution->push_back(value == Minisat::l_True ? 1 : 0);
		}
	}

	learnts = solver.get_learnts();

	return result;
}

std::vector<std::vector<std::pair<int, bool>>> MinisatInternalSatSolver::get_learnts() {
	return std::move(learnts);
}

bool MinisatInternalSatSolver::solve(const std::vector<Literal>& assumptions,
                                     std::vector<int>* solution) {
	Minisat::vec<Minisat::Lit> minisat_assumptions;
	for (auto& literal : assumptions) {
		if (literal.second) {
			minisat_assumptions.push(literals.at(literal.first));
		} else {
			minisat_assumptions.push(~literals.at(literal.first));
		}
	}
	bool result = solver.solve(minisat_assumptions, do_simplification);
	if (result) {
		solution->clear();
		for (auto& lit : literals) {
			auto value = solver.modelValue(lit);
			minimum_core_assert(value == Minisat::l_True || value == Minisat::l_False);
			solution->push_back(value == Minisat::l_True ? 1 : 0);
		}
	}
	return result;
}

void IP::save_CNF(const std::string& file_name) {
	auto minisat = new MinisatInternalSatSolver(true);
	std::unique_ptr<SatSolver> minisat_ptr(minisat);
	IpToSatSolver solver([&minisat_ptr]() { return std::move(minisat_ptr); });
	solver.convert(this);

	auto f = std::fopen(file_name.c_str(), "w");
	check(f != nullptr, "Could not open file.");
	at_scope_exit(std::fclose(f););
	minisat->get_solver().toDimacs(f, {});
}
}  // namespace linear
}  // namespace minimum
