
#include <algorithm>
#include <cmath>
using namespace std;

#include <minimum/algorithms/knapsack.h>
#include <minimum/core/check.h>
#include <iostream>
namespace minimum {
namespace algorithms {

template <typename Int>
Int gcd(Int a, Int b) {
	minimum_core_assert(b >= 0);
	if (a < b) {
		return gcd(b, a);
	}
	if (a == 1 || b == 1) {
		return 1;
	}
	if (b == 0) {
		return a;
	}
	return gcd(b, a - b);
}

double solve_knapsack(int weight_limit,
                      const vector<int>& weights,
                      const vector<double>& values,
                      vector<int>* solution) {
	int factor = weight_limit;
	for (auto w : weights) {
		factor = gcd(factor, w);
	}
	weight_limit /= factor;

	vector<double> m(weight_limit + 1);
	m[0] = 0;
	for (int w = 1; w <= weight_limit; ++w) {
		m[w] = 0;
		for (int i = 0; i < weights.size(); ++i) {
			auto weights_i = weights[i] / factor;
			if (weights_i <= w) {
				m[w] = max(m[w], values[i] + m[w - weights_i]);
			}
		}
	}

	solution->clear();
	solution->resize(weights.size(), 0);
	for (int w = weight_limit; w > 0;) {
		bool found = false;
		for (int i = 0; i < weights.size(); ++i) {
			auto weights_i = weights[i] / factor;
			if (weights_i <= w) {
				if (m[w] == values[i] + m[w - weights_i]) {
					solution->at(i) += 1;
					w -= weights_i;
					found = true;
					break;
				}
			}
		}
		if (!found) {
			--w;
		}
	}

	// Check solution. This is not needed but is very quick compared
	// to the code above.
	double solution_value = 0;
	int solution_weight = 0;
	for (int i = 0; i < weights.size(); ++i) {
		auto weights_i = weights[i] / factor;
		solution_value += solution->at(i) * values[i];
		solution_weight += solution->at(i) * weights_i;
	}
	minimum_core_assert(solution_weight <= weight_limit, solution_weight, " > ", weight_limit);
	minimum_core_assert(
	    abs(solution_value - m[weight_limit]) < 1e-4, solution_value, " != ", m[weight_limit]);

	return factor * m[weight_limit];
}
}  // namespace algorithms
}  // namespace minimum
