#include <catch.hpp>

#include <minimum/core/combinatorics.h>
using namespace minimum::core;
using namespace std;

bool contains(std::vector<std::vector<int>> subsets, std::vector<int> set) {
	std::sort(set.begin(), set.end());
	for (auto& subset : subsets) {
		std::sort(subset.begin(), subset.end());
	}
	return std::find(subsets.begin(), subsets.end(), set) != subsets.end();
}

TEST_CASE("subsets_size_1") {
	vector<int> s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	subsets.emplace_back(vector<int>({1, 2, 5}));  // To test that it is cleared.
	minimum::core::generate_subsets(s, 1, &subsets);
	CHECK(subsets.size() == 4);
	CHECK(contains(subsets, {1}));
	CHECK(contains(subsets, {2}));
	CHECK(contains(subsets, {5}));
	CHECK(contains(subsets, {8}));
}

TEST_CASE("subsets_size_2") {
	vector<int> s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	generate_subsets(s, 2, &subsets);
	CHECK(subsets.size() == 4 * 3 / 2);
	CHECK(contains(subsets, {1, 2}));
	CHECK(contains(subsets, {1, 5}));
	CHECK(contains(subsets, {1, 8}));
	CHECK(contains(subsets, {5, 2}));
	CHECK(contains(subsets, {2, 8}));
	CHECK(contains(subsets, {5, 8}));
}

TEST_CASE("subsets_size_3") {
	vector<int> s = {1, 5, 7, 8, 9};
	vector<vector<int>> subsets;
	generate_subsets(s, 3, &subsets);
	CHECK(subsets.size() == 10);
	for (const auto& subset : subsets) {
		CHECK(subset.size() == 3);
	}
}

TEST_CASE("subsets_size_n") {
	vector<int> s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	generate_subsets(s, 4, &subsets);
	CHECK(subsets.size() == 1);
	CHECK(contains(subsets, {1, 2, 5, 8}));
}

TEST_CASE("subsets_size_15_k_7") {
	vector<int> s = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	vector<vector<int>> subsets;
	generate_subsets(s, 7, &subsets);
	CHECK(subsets.size() == 6435);
}

TEST_CASE("too_many_subsets") {
	vector<int> s(50);
	vector<vector<int>> subsets;
	CHECK_THROWS(generate_subsets(s, 25, &subsets));
}
