#include <cmath>
#include <memory>
#include <vector>

#include "projection_solver.h"

#include <minimum/core/range.h>
#include <minimum/linear/projection_solver.h>
#include <minimum/linear/proto.h>
using namespace minimum::core;
using namespace std;

namespace minimum {
namespace linear {
class ProjectionSolutions : public Solutions {
   public:
	ProjectionSolutions(IP* ip_, const ProjectionSolver& solver_)
	    : ip(ip_), solver(solver_), m(ip_->get().constraint_size()), n(ip_->get().variable_size()) {
		variable_bound.reserve(n);
		for (auto j : range(n)) {
			auto& bound = ip->get().variable(j).bound();
			variable_bound.emplace_back(bound.lower(), bound.upper());
		}
		int start = 0;
		for (auto& constraint : ip->get().constraint()) {
			auto& bound = constraint.bound();
			if (abs(bound.lower() - bound.upper()) < 1e-9) {
				double norm_squared = 0;
				for (auto& entry : constraint.sum()) {
					equality_entries.emplace_back(entry.variable(), entry.coefficient());
					norm_squared += entry.coefficient() * entry.coefficient();
				}
				double norm = sqrt(norm_squared);
				for (int idx = start; idx < equality_entries.size(); ++idx) {
					equality_entries[idx].coefficient /= norm;
				}
				equality_info.emplace_back(start, bound.lower() / norm);
				start = equality_entries.size();
			}
		}
		equality_info.emplace_back(start, numeric_limits<double>::quiet_NaN());

		start = 0;
		for (auto& constraint : ip->get().constraint()) {
			auto& bound = constraint.bound();
			if (abs(bound.lower() - bound.upper()) >= 1e-9) {
				check(bound.lower() <= -1e100 || bound.upper() >= 1e100, "One bound needed.");
				double swap = 1.0;
				double rhs = bound.upper();
				if (bound.upper() >= 1e100) {
					swap = -1.0;
					rhs = bound.lower();
				}
				double norm_squared = 0;
				for (auto& entry : constraint.sum()) {
					inequality_entries.emplace_back(entry.variable(), entry.coefficient());
					norm_squared += entry.coefficient() * entry.coefficient();
				}
				double norm = sqrt(norm_squared);
				for (int idx = start; idx < inequality_entries.size(); ++idx) {
					inequality_entries[idx].coefficient /= swap * norm;
				}
				inequality_info.emplace_back(start, swap * rhs / norm);
				start = inequality_entries.size();
			}
		}
		inequality_info.emplace_back(start, numeric_limits<double>::quiet_NaN());
	}

	virtual bool get() override {
		vector<double> x(n, 0);

		for (int iter = 1; iter <= solver.max_iterations; ++iter) {
			for (auto j : range(n)) {
				x[j] = max(x[j], variable_bound[j].first);
				x[j] = min(x[j], variable_bound[j].second);
			}

			for (int i = 0; i < int(equality_info.size()) - 1; ++i) {
				auto rhs = equality_info[i].rhs;
				double aTx = 0;
				for (int idx = equality_info[i].start; idx < equality_info[i + 1].start; ++idx) {
					auto& entry = equality_entries[idx];
					aTx += entry.coefficient * x[entry.index];
				}
				auto diff = aTx - rhs;
				if (abs(diff) > 1e-9) {
					for (int idx = equality_info[i].start; idx < equality_info[i + 1].start;
					     ++idx) {
						auto& entry = equality_entries[idx];
						auto& xj = x[entry.index];
						xj = xj - diff * entry.coefficient;
					}
				}
			}

			for (int i = 0; i < int(inequality_info.size()) - 1; ++i) {
				auto rhs = inequality_info[i].rhs;
				double aTx = 0;
				for (int idx = inequality_info[i].start; idx < inequality_info[i + 1].start;
				     ++idx) {
					auto& entry = inequality_entries[idx];
					aTx += entry.coefficient * x[entry.index];
				}
				auto diff = aTx - rhs;
				if (diff > 0) {
					for (int idx = inequality_info[i].start; idx < inequality_info[i + 1].start;
					     ++idx) {
						auto& entry = inequality_entries[idx];
						auto& xj = x[entry.index];
						xj = xj - diff * entry.coefficient;
					}
				}
			}

			if (iter % solver.test_interval == 0) {
				clog << "Iteration=" << iter << "\n";
				for (auto j : range(n)) {
					ip->set_solution(j, x[j]);
				}
				if (ip->is_feasible(solver.tolerance)) {
					return true;
				}
			}
		}
		return false;
	}

   private:
	IP* ip;
	const ProjectionSolver& solver;
	int m;
	int n;
	vector<pair<double, double>> variable_bound;

	struct ConstraintInfo {
		ConstraintInfo(int start_, double rhs_) : start(start_), rhs(rhs_) {}
		int start;
		double rhs;
	};
	struct Entry {
		Entry(int index_, double coefficient_) : index(index_), coefficient(coefficient_) {}
		int index;
		double coefficient;
	};
	vector<ConstraintInfo> equality_info;
	vector<Entry> equality_entries;
	vector<ConstraintInfo> inequality_info;
	vector<Entry> inequality_entries;
};

SolutionsPointer ProjectionSolver::solutions(IP* ip) const {
	unique_ptr<Solutions> solutions = make_unique<ProjectionSolutions>(ip, *this);
	return move(solutions);
}
}  // namespace linear
}  // namespace minimum
