#pragma once

#include <vector>

namespace minimum {
namespace core {

// Creates a list of all subsets of a given list.
template <typename T>
void generate_subsets(const std::vector<T>& set,
                      int subset_size,
                      std::vector<std::vector<T>>* output);

namespace internal {
template <typename T>
void subset(const std::vector<T>& set,
            int left,
            int index,
            std::vector<T>* scratch_space,
            std::vector<std::vector<T>>* all_subsets) {
	if (left == 0) {
		all_subsets->push_back(*scratch_space);
		return;
	}
	if (left > set.size() - index) {
		// We don’t have enough elements left to create a subset.
		return;
	}
	for (std::size_t i = index; i < set.size(); i++) {
		scratch_space->push_back(set[i]);
		subset(set, left - 1, i + 1, scratch_space, all_subsets);
		scratch_space->pop_back();
	}
}

size_t choose(size_t n, size_t k) {
	if (k == 0) {
		return 1;
	}
	return (n * choose(n - 1, k - 1)) / k;
}
}  // namespace internal.

template <typename T>
void generate_subsets(const std::vector<T>& set,
                      int subset_size,
                      std::vector<std::vector<T>>* output) {
	size_t num_subsets = internal::choose(set.size(), subset_size);
	if (num_subsets > 50000000) {
		// Maybe change this limit in the future.
		throw std::runtime_error("Too many subsets. Choose a better algorithm.");
	}

	output->clear();
	output->reserve(num_subsets);
	std::vector<T> scratch_space;
	scratch_space.reserve(subset_size);
	internal::subset(set, subset_size, 0, &scratch_space, output);
}
}  // namespace core
}  // namespace minimum
