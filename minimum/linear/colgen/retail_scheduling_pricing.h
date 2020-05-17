
#include <random>
#include <vector>

#include <minimum/linear/colgen/export.h>
#include <minimum/linear/retail_scheduling.h>

namespace minimum {
namespace linear {
namespace colgen {

bool MINIMUM_LINEAR_COLGEN_API create_roster_graph(const RetailProblem& problem,
                                                   const std::vector<double>& duals,
                                                   int staff_index,
                                                   const std::vector<std::vector<int>>& fixes,
                                                   std::vector<std::vector<int>>* solution,
                                                   std::mt19937* rng);
}
}  // namespace linear
}  // namespace minimum
