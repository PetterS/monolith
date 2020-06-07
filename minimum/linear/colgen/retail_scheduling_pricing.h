
#include <functional>
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

struct RetailLocalSearchInfo {
	int computed_objective_value = -1;
	double elapsed_time = -1;
};

struct RetailLocalSearchParameters {
	// If the callback returns false, the search will stop.
	std::function<bool(const RetailLocalSearchInfo& info,
	                   const std::vector<std::vector<std::vector<int>>>& solution)>
	    callback = nullptr;
	// Set time limit to negative to disable time limit.
	double time_limit_seconds = 10.0;

	// If non-empty, used to initialize the search.
	std::vector<std::vector<std::vector<int>>> solution;
};

void MINIMUM_LINEAR_COLGEN_API retail_local_search(const RetailProblem& problem,
                                                   const RetailLocalSearchParameters& parameters);

// Available for testing.
std::vector<int> MINIMUM_LINEAR_COLGEN_API
create_retail_shift(const std::vector<std::vector<double>>& shift_costs);

}  // namespace colgen
}  // namespace linear
}  // namespace minimum
