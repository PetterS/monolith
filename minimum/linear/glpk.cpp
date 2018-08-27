
#include <cmath>
#include <unordered_map>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/scope_guard.h>
#include <minimum/linear/easy-ip-internal.h>
#include <minimum/linear/glpk.h>
#include <minimum/linear/ip.h>

using namespace minimum::core;

namespace std {
template <typename T, typename U>
struct hash<pair<T, U>> {
	size_t operator()(pair<T, U> p) const { return hash<T>()(p.first) ^ hash<U>()(p.second); }
};
}  // namespace std

namespace minimum {
namespace linear {

class GlpkSolutions : public Solutions {
   public:
	GlpkSolutions(IP* ip_, bool silent) : ip(ip_) {
		ip->linearize_pseudoboolean_terms();

		problem.add_columns(ip->get_number_of_variables());
		problem.add_rows(ip->get().constraint_size());
		problem.set_cost(*ip);
		problem.set_constant(ip->get_objective_constant());
		problem.set_col_bounds(*ip);
		problem.set_row_bounds(*ip);
		problem.load_matrix(*ip);

		for (auto j : range(ip->get().variable_size())) {
			if (ip->get().variable(j).type() == proto::Variable_Type_INTEGER) {
				problem.set_col_integer(j, true);
			}
		}

		if (silent) {
			problem.set_silent();
		}
	}

	virtual bool get() override {
		if (first_time) {
			first_time = false;
			return solve();
		} else {
			return next();
		}
	}

	bool solve() {
		bool success = problem.solve();
		if (success) {
			problem.extract_solution(ip);
		}
		return success;
	}

	bool next() {
		//
		// We add two new rows to the problem in order to get the
		// next solution. If the current solution is x = (1, 0, 1),
		//
		//    (a)  (1 - x1) + x2 + (1 - x3) >= 1
		//    (b)  objective(x) == *optimal*.
		//
		vector<int> constraint_row_indices;
		vector<double> constraint_row_values;
		constraint_row_indices.push_back(-1);
		constraint_row_values.push_back(-1);

		vector<int> objective_constraint_row_indices;
		vector<double> objective_constraint_row_values;
		objective_constraint_row_indices.push_back(-1);
		objective_constraint_row_values.push_back(-1);

		double best_objective = 0;
		double solution_rhs = 1.0;
		bool found_integer_variable = false;
		for (int iColumn = 0; iColumn < ip->get().variable_size(); iColumn++) {
			auto& variable = ip->get().variable(iColumn);
			double value = ip->get_solution().primal(iColumn);
			if (variable.type() == proto::Variable_Type_INTEGER) {
				found_integer_variable = true;
			}
			if (variable.type() == proto::Variable_Type_INTEGER && !variable.is_helper()) {
				// only works for 0-1 variables
				auto& bound = ip->get().variable(iColumn).bound();
				minimum_core_assert(bound.lower() == 0.0 || bound.lower() == 1.0);
				minimum_core_assert(bound.upper() == 0.0 || bound.upper() == 1.0);
				// double check integer
				minimum_core_assert(fabs(floor(value + 0.5) - value) < 1.0e-5);

				constraint_row_indices.push_back(iColumn + 1);
				if (value > 0.5) {
					// at 1.0
					constraint_row_values.push_back(-1.0);
					solution_rhs -= 1.0;
				} else {
					// at 0.0
					constraint_row_values.push_back(1.0);
				}
			}

			best_objective += value * ip->get().variable(iColumn).cost();
			objective_constraint_row_indices.push_back(iColumn + 1);
			objective_constraint_row_values.push_back(ip->get().variable(iColumn).cost());
		}
		check(found_integer_variable, "Need a boolean variable.");

		glp_add_rows(problem.get_problem(), 2);
		int constraint_row = glp_get_num_rows(problem.get_problem());
		int solution_row = constraint_row - 1;
		glp_set_mat_row(problem.get_problem(),
		                constraint_row,
		                constraint_row_values.size() - 1,
		                constraint_row_indices.data(),
		                constraint_row_values.data());
		glp_set_row_bnds(problem.get_problem(), constraint_row, GLP_LO, solution_rhs, 999);
		glp_set_mat_row(problem.get_problem(),
		                solution_row,
		                objective_constraint_row_values.size() - 1,
		                objective_constraint_row_indices.data(),
		                objective_constraint_row_values.data());
		glp_set_row_bnds(
		    problem.get_problem(), solution_row, GLP_FX, best_objective, best_objective);

		bool success = problem.solve();
		if (success) {
			problem.extract_solution(ip);
		}

		return success;
	}

	GlpProblem problem;
	IP* ip;
	bool first_time = true;
};

static int contraint_type(double lower, double upper) {
	bool has_lower = lower > -1e100;
	bool has_upper = upper < 1e100;
	if (lower >= upper) {
		return GLP_FX;
	} else if (has_lower && has_upper) {
		return GLP_DB;
	} else if (has_lower) {
		return GLP_LO;
	} else if (has_upper) {
		return GLP_UP;
	} else {
		return GLP_FR;
	}
}

thread_local int environment_ref_count = 0;

GlpProblem::GlpProblem() : has_integers(false) {
	environment_ref_count++;

	problem = glp_create_prob();
	glp_set_obj_dir(problem, GLP_MIN);

	glp_init_smcp(&smcp);
	smcp.presolve = GLP_ON;

	glp_init_iocp(&iocp);
	iocp.presolve = GLP_ON;

	iocp.fp_heur = GLP_ON;
	// iocp.ps_heur = GLP_ON;
	iocp.gmi_cuts = GLP_ON;
	iocp.mir_cuts = GLP_ON;
	iocp.cov_cuts = GLP_ON;
	// iocp.clq_cuts = GLP_ON;  // Crashes with GCC for some tests.
}

GlpProblem::~GlpProblem() {
	if (problem != nullptr) {
		glp_delete_prob(problem);
		environment_ref_count--;
		minimum_core_assert(environment_ref_count >= 0);
		if (environment_ref_count == 0) {
			glp_free_env();
		}
	}
}

void GlpProblem::operator=(GlpProblem&& rhs) {
	problem = rhs.problem;
	iocp = rhs.iocp;
	smcp = rhs.smcp;
	has_integers = rhs.has_integers;
	rhs.problem = nullptr;
}

void GlpProblem::add_columns(int m) { glp_add_cols(problem, m); }

int GlpProblem::num_cols() const { return glp_get_num_cols(problem); }

void GlpProblem::set_col_bounds(const IP& ip) {
	for (auto j : range(ip.get().variable_size())) {
		auto& bound = ip.get().variable(j).bound();
		glp_set_col_bnds(problem,
		                 j + 1,
		                 contraint_type(bound.lower(), bound.upper()),
		                 bound.lower(),
		                 bound.upper());
	}
}

void GlpProblem::set_col_integer(int j, bool is_integer) {
	if (is_integer) {
		has_integers = true;
	}
	int kind = is_integer ? GLP_IV : GLP_CV;
	glp_set_col_kind(problem, j + 1, kind);
}

void GlpProblem::set_cost(const IP& ip) {
	for (auto j : range(ip.get().variable_size())) {
		glp_set_obj_coef(problem, j + 1, ip.get().variable(j).cost());
	}
}

void GlpProblem::set_constant(double cost) { glp_set_obj_coef(problem, 0, cost); }

void GlpProblem::add_rows(int n) { glp_add_rows(problem, n); }

int GlpProblem::num_rows() const { return glp_get_num_rows(problem); }

void GlpProblem::set_row_bounds(const IP& ip) {
	for (auto i : range(ip.get().constraint_size())) {
		auto& bound = ip.get().constraint(i).bound();
		glp_set_row_bnds(problem,
		                 i + 1,
		                 contraint_type(bound.lower(), bound.upper()),
		                 bound.lower(),
		                 bound.upper());
	}
}

void GlpProblem::load_matrix(const IP& ip) {
	auto& ip_data = ip.get();
	auto m = ip_data.constraint_size();
	auto matrix_size = ip.matrix_size();

	vector<int> glpk_rows, glpk_cols;
	vector<double> glpk_values;
	glpk_rows.reserve(matrix_size + 1);
	glpk_cols.reserve(matrix_size + 1);
	glpk_values.reserve(matrix_size + 1);
	glpk_rows.push_back(-1);
	glpk_cols.push_back(-1);
	glpk_values.push_back(-1);

	for (auto i : range(m)) {
		auto& constraint = ip_data.constraint(i);
		for (auto& entry : constraint.sum()) {
			glpk_values.push_back(entry.coefficient());
			glpk_cols.push_back(entry.variable() + 1);
			glpk_rows.push_back(i + 1);
		}
	}

	glp_load_matrix(
	    problem, glpk_rows.size() - 1, glpk_rows.data(), glpk_cols.data(), glpk_values.data());
}

bool GlpProblem::solve() {
	if (!has_integers) {
		// If we don’t have integer variables, just use simplex,
		// because we will then be able to extract the duals as
		// well.
		auto result = glp_simplex(problem, &smcp);
		if (result != 0) {
			return false;
		}
		return glp_get_status(problem) == GLP_OPT;
	}
	bool success = glp_intopt(problem, &iocp) == 0;
	return success && glp_mip_status(problem) == GLP_OPT;
}

void GlpProblem::set_silent() {
	iocp.msg_lev = GLP_MSG_OFF;
	smcp.msg_lev = GLP_MSG_OFF;
}

void GlpProblem::extract_solution(IP* ip) {
	if (!has_integers) {
		int m = glp_get_num_cols(problem);
		for (int j = 0; j < m; ++j) {
			ip->set_solution(j, glp_get_col_prim(problem, j + 1));
		}

		int n = glp_get_num_rows(problem);
		for (int i = 0; i < n; ++i) {
			ip->set_dual_solution(i, glp_get_row_dual(problem, i + 1));
		}
	} else {
		int m = glp_get_num_cols(problem);
		for (int j = 0; j < m; ++j) {
			ip->set_solution(j, glp_mip_col_val(problem, j + 1));
		}

		// TODO: Duals?
	}
}

void GlpProblem::set_presolve(bool enable) {
	if (enable) {
		smcp.presolve = GLP_ON;
	} else {
		smcp.presolve = GLP_OFF;
	}
}

GlpkSolver::GlpkSolver() : Solver() {}

SolutionsPointer GlpkSolver::solutions(IP* ip) const {
	return {std::make_unique<GlpkSolutions>(ip, silent)};
}
}  // namespace linear
}  // namespace minimum
