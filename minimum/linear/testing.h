#pragma once

#include <set>
#include <string>

#include <minimum/linear/export.h>
#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {

// Creates a Sudoku IP with no given numbers.
MINIMUM_LINEAR_API auto create_soduku_IP(IP& ip, int n = 3)
    -> decltype(ip.add_boolean_cube(9, 9, 9));

// Performs simple smoke test for linear programming.
MINIMUM_LINEAR_API void simple_linear_programming_tests(const Solver& solver,
                                                        std::set<std::string> tests_to_skip = {});

// Performs simple smoke tests for integer programming.
MINIMUM_LINEAR_API void simple_integer_programming_tests(const Solver& solver,
                                                         std::set<std::string> tests_to_skip = {});
}  // namespace linear
}  // namespace minimum
