#pragma once
#include <iostream>
#include <vector>

#include <minimum/linear/export.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/proto.h>

namespace minimum {
namespace linear {

// Verifies that a SchedulingProblem is valid. Throws if any error is discovered.
MINIMUM_LINEAR_API void verify_scheduling_problem(
    const minimum::linear::proto::SchedulingProblem& problem);

MINIMUM_LINEAR_API minimum::linear::proto::SchedulingProblem read_nottingham_instance(
    std::istream& in, bool verbose = true);

MINIMUM_LINEAR_API minimum::linear::proto::SchedulingProblem read_gent_instance(
    std::istream& problem_data, std::istream& case_data);

MINIMUM_LINEAR_API minimum::linear::proto::TaskSchedulingProblem read_ptask_instance(
    std::istream& data);

MINIMUM_LINEAR_API void save_solution(std::ostream& out,
                                      const minimum::linear::proto::SchedulingProblem& problem,
                                      const std::vector<std::vector<std::vector<int>>>& solution);

MINIMUM_LINEAR_API std::vector<std::vector<std::vector<int>>> load_solution(
    std::istream& in, const minimum::linear::proto::SchedulingProblem& problem);

MINIMUM_LINEAR_API double print_solution(
    std::ostream& out,
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<std::vector<std::vector<int>>>& solution);

MINIMUM_LINEAR_API int objective_value(
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<std::vector<std::vector<int>>>& current_solution);

MINIMUM_LINEAR_API int cover_cost(const minimum::linear::proto::SchedulingProblem& problem,
                                  const std::vector<std::vector<std::vector<int>>>& solution);

MINIMUM_LINEAR_API minimum::linear::Sum create_ip(
    minimum::linear::IP& ip,
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<std::vector<std::vector<minimum::linear::BooleanVariable>>>& x);

MINIMUM_LINEAR_API minimum::linear::Sum create_ip(
    minimum::linear::IP& ip,
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<std::vector<std::vector<minimum::linear::Sum>>>& x);

MINIMUM_LINEAR_API minimum::linear::Sum create_ip(
    minimum::linear::IP& ip,
    const minimum::linear::proto::SchedulingProblem& problem,
    const std::vector<std::vector<std::vector<minimum::linear::Sum>>>& x,
    const std::vector<int>& staff_indices,
    const std::vector<int>& day_indices);

// Tries to improve the solution by simply assigning more shifts to personnel.
// Returns a description of the operation if successdul, otherwise an empty string.
MINIMUM_LINEAR_API std::string quick_solution_improvement(
    const proto::SchedulingProblem& problem, vector<vector<vector<int>>>& current_solution);

}  // namespace linear
}  // namespace minimum