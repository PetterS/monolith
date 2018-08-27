#pragma once
#include <glpk.h>

#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

// GlpProblem is not safe to use from multiple threads at the same time.
// This is a limitation of Glpk.
class MINIMUM_LINEAR_API GlpProblem {
   public:
	GlpProblem();
	~GlpProblem();
	void operator=(GlpProblem&& rhs);

	void add_columns(int m);
	int num_cols() const;
	void set_col_bounds(const IP& ip);
	void set_col_integer(int j, bool is_integer);
	void set_cost(const IP& ip);
	void set_constant(double cost);

	void add_rows(int n);
	int num_rows() const;
	void set_row_bounds(const IP& ip);

	void load_matrix(const IP& ip);

	glp_prob* get_problem() { return problem; }

	bool solve();
	void set_silent();
	void extract_solution(IP* ip);

	// If presolve is disabled, the solver can be warmstarted.
	// Default is enabled.
	void set_presolve(bool enable);

   private:
	glp_prob* problem;
	glp_iocp iocp;
	glp_smcp smcp;
	bool has_integers;
};

class MINIMUM_LINEAR_API GlpkSolver : public Solver {
   public:
	GlpkSolver();
	virtual SolutionsPointer solutions(IP* ip) const override;
};
}  // namespace linear
}  // namespace minimum
