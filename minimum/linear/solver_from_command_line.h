#pragma once

#include <memory>

#include <minimum/linear/export.h>
#include <minimum/linear/solver.h>

namespace minimum {
namespace linear {
// Creates a solver for a given IP. The solver used is determined via
// the command line option --solver.
MINIMUM_LINEAR_API std::unique_ptr<Solver> create_solver_from_command_line();
}  // namespace linear
}  // namespace minimum
