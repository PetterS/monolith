#include <third-party/cleaneling/cleaneling.h>

#include <minimum/core/check.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/sat_solver.h>

using minimum::core::check;

namespace minimum {
namespace linear {

class CleanelingSolver : public SatSolver {
   public:
	CleanelingSolver() { solver = Cleaneling::new_solver(); }
	~CleanelingSolver() { delete solver; }

	virtual Var add_variable() override;
	virtual void add_clause(const std::vector<Literal>& clause) override;
	virtual bool solve(std::vector<int>* solution) override;
	virtual bool solve(const std::vector<Literal>& assumptions,
	                   std::vector<int>* solution) override;

   private:
	// SimpSolver is faster on Nurikabe, but hangs the unit tests with output.
	Cleaneling::AbstractSolver* solver;
	int var_count = 0;
};

std::unique_ptr<SatSolver> cleaneling_solver() { return std::make_unique<CleanelingSolver>(); }

SatSolver::Var CleanelingSolver::add_variable() {
	var_count++;
	return var_count - 1;
}

void CleanelingSolver::add_clause(const std::vector<Literal>& clause) {
	for (auto& literal : clause) {
		if (literal.second) {
			solver->addlit(literal.first + 1);
		} else {
			solver->addlit(-(literal.first + 1));
		}
	}
	solver->addlit(0);
}

bool CleanelingSolver::solve(std::vector<int>* solution) {
	auto result = solver->solve();
	if (result == Cleaneling::SATISFIABLE) {
		solution->clear();
		for (int i = 0; i < var_count; ++i) {
			auto value = solver->value(i + 1);
			// I don’t know why UNASSIGNED would be returned here.
			// minimum_core_assert(value == Cleaneling::TRUE || value == Cleaneling::FALSE);
			solution->push_back(value == Cleaneling::TRUE ? 1 : 0);
		}
	}
	return result == Cleaneling::SATISFIABLE;
}

bool CleanelingSolver::solve(const std::vector<Literal>& assumptions, std::vector<int>* solution) {
	check(false, "Solving with assumptions is not supported for Cleaneling.");
	return false;
}
}  // namespace linear
}  // namespace minimum
