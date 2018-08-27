#pragma once
#include <vector>

#include <minimum/isomorphism/graph-isomorphism.h>
#include <minimum/linear/ip.h>

namespace minimum {
namespace linear {

struct Automorphism {
	std::vector<minimum::isomorphism::Permutation> generators;
	std::vector<std::vector<int>> orbits;
};

MINIMUM_LINEAR_API Automorphism analyze_symmetry(const IP& ip);

// Experimental test code.
MINIMUM_LINEAR_API void remove_symmetries(IP* ip);
}  // namespace linear
}  // namespace minimum
