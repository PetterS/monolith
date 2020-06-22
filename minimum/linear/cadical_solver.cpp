#include <iterator>
#include <string>

#include <third-party/cadical/cadical/cadical.hpp>

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

class CadicalInternalSatSolver : public SatSolver {
   public:
	CadicalInternalSatSolver(bool silent_mode);

	virtual Var add_variable() override;
	virtual Var add_helper_variable() override;
	virtual void add_clause(const std::vector<Literal>& clause) override;
	virtual bool solve(std::vector<int>* solution) override;
	virtual bool solve(const std::vector<Literal>& assumptions,
	                   std::vector<int>* solution) override;

   private:
	CaDiCaL::Solver solver;
	std::vector<int> literals;
};

std::unique_ptr<SatSolver> cadical_solver(bool silent_mode) {
	return std::make_unique<CadicalInternalSatSolver>(silent_mode);
}

CadicalInternalSatSolver::CadicalInternalSatSolver(bool silent_mode) {}

SatSolver::Var CadicalInternalSatSolver::add_variable() {
	literals.push_back(literals.size() + 1);
	return literals.back() - 1;
}

SatSolver::Var CadicalInternalSatSolver::add_helper_variable() { return add_variable(); }

void CadicalInternalSatSolver::add_clause(const std::vector<Literal>& clause) {
	for (auto& literal : clause) {
		if (literal.second) {
			solver.add(literal.first + 1);
		} else {
			solver.add(-literal.first - 1);
		}
	}
	solver.add(0);
}

bool CadicalInternalSatSolver::solve(std::vector<int>* solution) {
	bool result = solver.solve() == 10;
	if (result) {
		solution->clear();
		for (auto& lit : literals) {
			auto value = solver.val(lit);
			solution->push_back(value > 0 ? 1 : 0);
		}
	}
	return result;
}

bool CadicalInternalSatSolver::solve(const std::vector<Literal>& assumptions,
                                     std::vector<int>* solution) {
	for (auto& literal : assumptions) {
		if (literal.second) {
			solver.assume(literal.first + 1);
		} else {
			solver.assume(-literal.first - 1);
		}
	}
	return solve(solution);
}
}  // namespace linear
}  // namespace minimum
