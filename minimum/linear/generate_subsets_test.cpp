// Petter Strandmark 2014.
#include <algorithm>
#include <catch.hpp>
#include <iostream>

#include <minimum/linear/ip.h>

using minimum::linear::generate_subsets;

bool contains(std::vector<std::vector<int>> subsets, std::vector<int> set) {
	std::sort(set.begin(), set.end());
	for (auto& subset : subsets) {
		std::sort(subset.begin(), subset.end());
	}
	return std::find(subsets.begin(), subsets.end(), set) != subsets.end();
}

void print(std::vector<std::vector<int>> subsets) {
	using namespace std;
	clog << "{";
	for (auto set : subsets) {
		clog << "{";
		for (auto elem : set) {
			clog << elem << " ";
		}
		clog << "} ";
	}
	clog << "}" << endl;
}

using namespace std;

TEST_CASE("subsets_size_1") {
	auto s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	subsets.emplace_back(vector<int>({1, 2, 5}));  // To test that it is cleared.
	generate_subsets(s, 1, &subsets);
	print(subsets);
	CHECK(subsets.size() == 4);
	CHECK(contains(subsets, {1}));
	CHECK(contains(subsets, {2}));
	CHECK(contains(subsets, {5}));
	CHECK(contains(subsets, {8}));
}

TEST_CASE("subsets_size_2") {
	auto s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	generate_subsets(s, 2, &subsets);
	print(subsets);
	CHECK(subsets.size() == 4 * 3 / 2);
	CHECK(contains(subsets, {1, 2}));
	CHECK(contains(subsets, {1, 5}));
	CHECK(contains(subsets, {1, 8}));
	CHECK(contains(subsets, {5, 2}));
	CHECK(contains(subsets, {2, 8}));
	CHECK(contains(subsets, {5, 8}));
}

TEST_CASE("subsets_size_3") {
	auto s = {1, 5, 7, 8, 9};
	vector<vector<int>> subsets;
	generate_subsets(s, 3, &subsets);
	print(subsets);
	CHECK(subsets.size() == 10);
	for (const auto& subset : subsets) {
		CHECK(subset.size() == 3);
	}
}

TEST_CASE("subsets_size_n") {
	auto s = {1, 2, 5, 8};
	vector<vector<int>> subsets;
	generate_subsets(s, 4, &subsets);
	print(subsets);
	CHECK(subsets.size() == 1);
	CHECK(contains(subsets, {1, 2, 5, 8}));
}

TEST_CASE("subsets_size_15_k_7") {
	auto s = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	vector<vector<int>> subsets;
	generate_subsets(s, 7, &subsets);
	CHECK(subsets.size() == 6435);
}
