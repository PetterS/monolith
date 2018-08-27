// http://homepages.cae.wisc.edu/~linderot/talks/silo.pdf
// http://wpweb2.tepper.cmu.edu/fmargot/PDF/grpsurv.pdf

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/isomorphism/graph-isomorphism.h>
#include <minimum/linear/symmetry.h>
using namespace minimum::core;
using namespace minimum::isomorphism;

vector<vector<int>> compute_orbits(int n, const vector<Permutation>& generators) {
	// See e.g.
	// http://www.math.cornell.edu/~kbrown/6310/computation.pdf

	vector<vector<int>> orbits;
	vector<bool> done(n, false);

	for (int i = 0; i < n; ++i) {
		if (done[i]) {
			// Orbit already computed.
			continue;
		}

		// Compute the orbit for i.
		orbits.emplace_back();
		auto& orbit = orbits.back();
		orbit.push_back(i);

		for (int k = 0; k < orbit.size(); ++k) {
			int y = orbit[k];
			for (auto& permutation : generators) {
				int perm_y = permutation.permute(y);
				auto itr = find(orbit.begin(), orbit.end(), perm_y);
				if (itr == orbit.end()) {
					orbit.push_back(perm_y);
				}
			}
		}

		// Set the orbit for all elements in the orbit.
		for (auto y : orbit) {
			done[y] = true;
		}

		std::sort(orbit.begin(), orbit.end());
	}

	return orbits;
}

namespace minimum {
namespace linear {

Automorphism analyze_symmetry(const IP& ip) {
	bool verbose = false;

	auto n = ip.get_number_of_variables();
	auto m = ip.get().constraint_size();

	Eigen::SparseMatrix<size_t> graph(m + n, m + n);

	using VariableClass = tuple<double, double>;
	// Maps every variable type to a unique, non-negative integer.
	map<VariableClass, int> variable_classes;

	for (int j = 0; j < n; ++j) {
		auto& bound = ip.get().variable(j).bound();
		double lb = bound.lower();
		double ub = bound.upper();
		// TODO: Use variable type.

		auto sz = variable_classes.size();
		auto res = variable_classes.emplace(make_tuple(lb, ub), sz);
		auto type_index = res.first->second;

		for (int c = 0; c < type_index; ++c) {
			add_edge(graph, j, j);
		}
	}

	// Maps every coefficient to a unique, non-negative integer.
	map<double, int> coefficient_classes;
	vector<vector<pair<int, int>>> constraints(m);
	for (auto i : range(m)) {
		auto& constraint = ip.get().constraint(i);
		for (auto& entry : constraint.sum()) {
			auto sz = coefficient_classes.size();
			auto res = coefficient_classes.emplace(entry.coefficient(), sz + 1);
			auto value_index = res.first->second;

			constraints[i].emplace_back(value_index, entry.variable());
		}
	}

	using ConstraintClass = tuple<double, double>;
	map<ConstraintClass, int> constraint_classes;
	for (int i = 0; i < m; ++i) {
		double lb = ip.get().constraint(i).bound().lower();
		double ub = ip.get().constraint(i).bound().lower();

		auto sz = constraint_classes.size();
		auto res = constraint_classes.emplace(make_tuple(lb, ub), sz);
		auto type_index = res.first->second;

		auto index = n + i;

		for (int c = 0; c < type_index; ++c) {
			add_edge(graph, index, index);
		}

		for (auto& entry : constraints[i]) {
			for (int c = 0; c < entry.first; ++c) {
				add_edge(graph, entry.second, index);
			}
		}
	}

	// TODO: Objective function.

	if (verbose) {
		cerr << "Graph has " << n + m << " nodes.\n";
		cerr << variable_classes.size() << " type(s) of variables: " << to_string(variable_classes)
		     << "\n";
		cerr << coefficient_classes.size()
		     << " type(s) of coefficients: " << to_string(coefficient_classes) << "\n";
		cerr << constraint_classes.size()
		     << " type(s) of constraints: " << to_string(constraint_classes) << "\n";
	}

	auto generators = bliss_automorphism(graph);
	auto orbits = compute_orbits(n, generators);

	if (verbose) {
		cerr << generators.size() << " generators. Orbits:\n";
		for (auto& orbit : orbits) {
			if (orbit.size() > 1) {
				cerr << to_string(orbit) << "\n";
			}
		}
	}

	Automorphism automorphism;
	automorphism.generators = move(generators);
	automorphism.orbits = move(orbits);

	return automorphism;
}

void remove_symmetries(IP* ip) {
	// Experimental test code.

	auto automorphism = analyze_symmetry(*ip);
	int n = ip->get_number_of_variables();
	for (auto& generator : automorphism.generators) {
		Sum p_sum;
		Sum m_sum;
		for (int i = 0; i < n; ++i) {
			auto xi = ip->get_variable(i);
			auto yi = ip->get_variable(generator.permute(i));
			if (i == 0) {
				ip->add_constraint(xi <= yi);
			} else {
				auto xim1 = ip->get_variable(i - 1);
				auto yim1 = ip->get_variable(generator.permute(i - 1));
				auto p = ip->add_boolean();
				auto m = ip->add_boolean();
				ip->add_constraint(p + m <= 1);
				ip->add_constraint(xim1 - yim1 == p - m);
				// p + m = 0 ⇔ xim1 = yim1
				p_sum += p;
				m_sum += m;
				ip->add_constraint(xi <= yi + p_sum + m_sum);
			}
		}
	}
}
}  // namespace linear
}  // namespace minimum
