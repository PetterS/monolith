#pragma once
#include <random>
#include <string>
#include <vector>

#include <minimum/linear/colgen/export.h>
#include <minimum/linear/proto.h>

namespace minimum {
namespace linear {
namespace colgen {
MINIMUM_LINEAR_COLGEN_API bool create_roster_cspp(
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<double>& dual_variables,
    int p,
    const std::vector<std::vector<int>>& fixes,
    std::vector<std::vector<int>>* solution_for_staff,
    std::mt19937_64* rng,
    const std::string& graph_output_file_name = "");
}
}  // namespace linear
}  // namespace minimum
