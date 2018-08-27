
#include <cmath>
#include <unordered_map>

// clang-format off
#include <scs/scs.h>
#include <scs/amatrix.h>
// clang-format on

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/scs.h>
using namespace minimum::core;
using namespace std;

namespace minimum {
namespace linear {

class ScsSolutions : public Solutions {
   public:
	ScsSolutions(IP* ip_, bool silent_);
	~ScsSolutions();

	virtual bool get() override;

	virtual std::string log() override { return solver_log; }

	void create_work() { work = scs_init(data, cone, &info); }

	void set_max_iterations(int max) { data->stgs->max_iters = max; }
	void set_convergence_tolerance(double eps) { data->stgs->eps = eps; }

	Data* data;
	Cone* cone;
	Sol* sol;
	Sol* opt_sol;
	Info info = {0};
	Work* work = nullptr;

	vector<int> constraint_indices;

	std::string solver_log;
	IP* ip;
	bool silent;
};

ScsSolutions::ScsSolutions(IP* ip_, bool silent_) : ip(ip_), silent(silent_) {
	cone = (Cone*)scs_calloc(1, sizeof(Cone));
	data = (Data*)scs_calloc(1, sizeof(Data));
	data->stgs = (Settings*)scs_calloc(1, sizeof(Settings));
	sol = (Sol*)scs_calloc(1, sizeof(Sol));
	opt_sol = (Sol*)scs_calloc(1, sizeof(Sol));

	setDefaultSettings(data);

	auto m = ip->get().constraint_size();
	auto n = ip->get_number_of_variables();
	data->n = n;
	data->c = (double*)scs_calloc(n, sizeof(double));

	for (int j = 0; j < n; j++) {
		auto& variable = ip->get().variable(j);
		auto& bound = variable.bound();
		if (bound.lower() > -1e100) {
			m++;
		}
		if (bound.upper() < 1e100) {
			m++;
		}
	}

	data->m = m;
	data->b = (double*)scs_calloc(m, sizeof(double));

	// Create permutation to put equality constraints first.
	vector<int> equality_indices;
	vector<int> inequality_indices;
	for (int i = 0; i < ip->get().constraint_size(); ++i) {
		auto& constraint = ip->get().constraint(i);
		auto& bound = constraint.bound();
		if (bound.lower() == bound.upper()) {
			equality_indices.push_back(i);
		} else {
			inequality_indices.push_back(i);
		}
	}
	constraint_indices.clear();
	constraint_indices.insert(
	    constraint_indices.end(), equality_indices.begin(), equality_indices.end());
	constraint_indices.insert(
	    constraint_indices.end(), inequality_indices.begin(), inequality_indices.end());

	// Create A and b matrix.
	int nnz = 0;
	cone->f = equality_indices.size();
	cone->l = m - equality_indices.size();
	vector<vector<pair<int, double>>> column_indices(n);
	for (int i = 0; i < ip->get().constraint_size(); ++i) {
		int org_row = constraint_indices[i];
		auto& constraint = ip->get().constraint(org_row);
		auto& bound = constraint.bound();

		double neg = 1.0;
		if (bound.upper() < 1e100) {
			// = or ≤ constraint.
			data->b[i] = bound.upper();
		} else {
			// ≥
			data->b[i] = -bound.lower();
			neg = -1.0;
		}

		for (auto& entry : constraint.sum()) {
			int j = entry.variable();
			column_indices[j].emplace_back(i, neg * entry.coefficient());
			nnz++;
		}
	}

	int next_constraint = ip->get().constraint_size();
	for (int j = 0; j < n; j++) {
		auto& variable = ip->get().variable(j);
		auto& bound = variable.bound();
		if (bound.upper() < 1e100) {
			column_indices[j].emplace_back(next_constraint, 1.0);
			data->b[next_constraint] = bound.upper();
			next_constraint++;
			nnz++;
		}
		if (bound.lower() > -1e100) {
			column_indices[j].emplace_back(next_constraint, -1.0);
			data->b[next_constraint] = -bound.lower();
			next_constraint++;
			nnz++;
		}

		data->c[j] = variable.cost();
	}
	minimum_core_assert(next_constraint == m);

	AMatrix* A = (AMatrix*)scs_calloc(1, sizeof(AMatrix));
	data->A = A;
	A->i = (int*)scs_calloc(nnz, sizeof(int));
	A->p = (int*)scs_calloc((n + 1), sizeof(int));
	A->x = (double*)scs_calloc(nnz, sizeof(double));
	A->n = n;
	A->m = m;
	int ind = 0;
	A->p[0] = 0;
	for (int j = 0; j < n; j++) {
		for (auto& entry : column_indices[j]) {
			auto i = entry.first; /* row */
			auto value = entry.second;
			A->x[ind] = value;
			A->i[ind] = i;
			ind++;
		}
		A->p[j + 1] = ind;
	}
	minimum_core_assert(ind == nnz, ind, " != ", nnz);

	sol->x = (double*)scs_calloc(data->n, sizeof(double));
	sol->y = (double*)scs_calloc(data->m, sizeof(double));
	sol->s = (double*)scs_calloc(data->m, sizeof(double));
}

ScsSolutions::~ScsSolutions() {
	if (work != nullptr) {
		scs_finish(work);
	}
	freeData(data, cone);
	freeSol(sol);
	freeSol(opt_sol);
}

ScsSolver::~ScsSolver() {}

bool ScsSolutions::get() {
	scs_also_print_to_stdout = silent ? 0 : 1;

	// TODO: Warm-start the solver with an initial guess. This is not
	// entirely trivial since we have added extra constraints for which
	// we do not have dual variables.
	if (!work) {
		create_work();
	}

	data->stgs->warm_start = 1;
	for (int i = 0; i < ip->get().constraint_size(); ++i) {
		int org_row = constraint_indices[i];
		auto dual = ip->get_solution().dual(org_row);
		if (dual != dual) {
			data->stgs->warm_start = 0;
		}
		if (ip->get().constraint(org_row).bound().lower()
		    == ip->get().constraint(org_row).bound().upper()) {
			// Negative sign for the duals in order to match the linear programming
			// duals.
			sol->y[i] = -dual;
		} else {
			sol->y[i] = dual;
		}

		// TODO: Warm-start s.
		sol->s[i] = 0;
	}
	// TODO: We have more y variables corresponding to box constraints.

	for (int j = 0; j < data->n; ++j) {
		auto value = ip->get_solution().primal(j);
		if (value != value) {
			data->stgs->warm_start = 0;
		}
		sol->x[j] = value;
	}

	scs_solve(work, data, cone, sol, &info);
	auto code = info.statusVal;

	for (int j = 0; j < data->n; ++j) {
		auto value = sol->x[j];
		auto& variable = ip->get().variable(j);

		// Make the bounds constraints always feasible.
		auto& bound = variable.bound();
		if (bound.upper() < 1e100) {
			value = min(value, bound.upper());
		}
		if (bound.lower() > -1e100) {
			value = max(value, bound.lower());
		}

		ip->set_solution(j, value);
	}

	for (int i = 0; i < ip->get().constraint_size(); ++i) {
		int org_row = constraint_indices[i];
		if (ip->get().constraint(org_row).bound().lower()
		    == ip->get().constraint(org_row).bound().upper()) {
			// Negative sign for the duals in order to match the linear programming
			// duals.
			ip->set_dual_solution(org_row, -sol->y[i]);
		} else {
			ip->set_dual_solution(org_row, sol->y[i]);
		}
	}

	solver_log = get_scs_log();
	return code == SCS_SOLVED || code == SCS_SOLVED_INACCURATE;
}

void ScsSolver::set_max_iterations(int max) { max_iterations = max; }

void ScsSolver::set_convergence_tolerance(double eps) { convergence_tolerance = eps; }

SolutionsPointer ScsSolver::solutions(IP* ip) const {
	auto sol = new ScsSolutions(ip, silent);
	if (max_iterations >= 0) {
		sol->set_max_iterations(max_iterations);
	}
	if (convergence_tolerance >= 0) {
		sol->set_convergence_tolerance(convergence_tolerance);
	}
	return {std::unique_ptr<Solutions>(sol)};
}
}  // namespace linear
}  // namespace minimum
