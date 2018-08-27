#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <vector>

#include <catch.hpp>

#include <minimum/core/immutable.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>

using namespace std;
using namespace minimum::core;

TEST_CASE("ImmutableSet_basic") {
	vector<int> v = {1, 5, 7, 9, 11};
	ImmutableSet<int> immutable_set(v.begin(), v.size());

	for (auto num : v) {
		CHECK(*immutable_set.find(num) == num);
		CHECK(immutable_set.contains(num));
		CHECK(immutable_set.count(num) == 1);
	}
	int count = 0;
	for (auto num : immutable_set) {
		count++;
	}
	CHECK(count == v.size());

	CHECK(immutable_set.find(-1) == immutable_set.end());
	CHECK(immutable_set.find(100) == immutable_set.end());
	CHECK(immutable_set.find(2) == immutable_set.end());
	CHECK(immutable_set.find(3) == immutable_set.end());
	CHECK(immutable_set.find(4) == immutable_set.end());
	CHECK(immutable_set.find(10) == immutable_set.end());
	CHECK(!immutable_set.contains(10));
	CHECK(immutable_set.count(10) == 0);
}

TEST_CASE("ImmutableSet_empty") {
	vector<int> v;
	ImmutableSet<int> immutable_set(v.begin(), 0);
	CHECK(immutable_set.find(10) == immutable_set.end());
}

TEST_CASE("ImmutableSet_exhaustive") {
	mt19937_64 rng;
	uniform_int_distribution<int> element_dist(0, 99);
	for (int size : range(1, 100)) {
		std::set<int> expected_set;
		for (int i : range(size)) {
			expected_set.insert(element_dist(rng));
		}

		ImmutableSet<int> immutable_set(expected_set.begin(), expected_set.size());
		for (int elem : range(-1, 101)) {
			CHECK((expected_set.find(elem) == expected_set.end())
			      == (immutable_set.find(elem) == immutable_set.end()));
		}
	}
}

TEST_CASE("ImmutableMap_basic") {
	map<int, string> m = {{1, "I"}, {5, "V"}, {7, "VII"}, {9, "IX"}, {11, "XI"}};
	ImmutableMap<int, string> immutable_map(m.begin(), m.size());

	REQUIRE(m.size() == immutable_map.size());
	for (auto entry : m) {
		CHECK((*immutable_map.find(entry.first)).second == m.at(entry.first));
		CHECK(immutable_map[entry.first] == entry.second);
		CHECK(immutable_map.at(entry.first) == entry.second);
		CHECK(immutable_map.count(entry.first) == 1);
	}
	int count = 0;
	for (auto entry : immutable_map) {
		count++;
	}
	CHECK(count == m.size());

	CHECK(immutable_map.find(-1) == immutable_map.end());
	CHECK(immutable_map.find(100) == immutable_map.end());
	CHECK(immutable_map.find(2) == immutable_map.end());
	CHECK(immutable_map.find(3) == immutable_map.end());
	CHECK(immutable_map.find(4) == immutable_map.end());
	CHECK(immutable_map.find(10) == immutable_map.end());
	CHECK_THROWS(immutable_map[-1]);
	CHECK_THROWS(immutable_map.at(-1));
	CHECK(immutable_map.count(-1) == 0);
}

TEST_CASE("ImmutableMap_iterator_arrow") {
	map<int, string> m = {{1, "I"}};
	ImmutableMap<int, string> immutable_map(m.begin(), m.size());
	auto itr = immutable_map.find(1);

	CHECK(itr->first == 1);
	CHECK(itr->second == "I");
}
