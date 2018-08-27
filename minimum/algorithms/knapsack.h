#pragma once
#include <vector>

#include <minimum/algorithms/export.h>

namespace minimum {
namespace algorithms {

// Simple knapsack solver that finds the highest value in O(weight_limit · weights.size()) time.
MINIMUM_ALGORITHMS_API double solve_knapsack(int weight_limit,
                                             const vector<int>& weights,
                                             const vector<double>& values,
                                             vector<int>* solution);
}  // namespace algorithms
}  // namespace minimum
