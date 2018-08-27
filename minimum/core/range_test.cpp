#include <array>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include <catch.hpp>

#include <minimum/core/range.h>
using namespace minimum::core;

TEST_CASE("single") {
	std::vector<int> v;
	for (int i : range(5)) {
		v.push_back(i);
	}
	CHECK(v == std::vector<int>({0, 1, 2, 3, 4}));
}

TEST_CASE("double") {
	std::vector<int> v;
	for (int i : range(4, 8)) {
		v.push_back(i);
	}
	CHECK(v == std::vector<int>({4, 5, 6, 7}));
}

TEST_CASE("negative") {
	std::vector<int> v;
	for (int i : range(-2, 1)) {
		v.push_back(i);
	}
	CHECK(v == std::vector<int>({-2, -1, 0}));
}

TEST_CASE("zip") {
	const vector<int> v1 = {1, 2, 3};
	vector<string> v2 = {"I", "II", "III"};

	vector<pair<int, string>> result;
	for (auto elem : zip(v1, v2)) {
		result.emplace_back(get<0>(elem), get<1>(elem));
	}
	vector<pair<int, string>> expected = {{1, "I"}, {2, "II"}, {3, "III"}};
	CHECK(result == expected);
}

template <typename T>
struct NoCopy {
	T t;
	NoCopy(T t_) : t(t_) {}
	NoCopy(const NoCopy&) = delete;
	NoCopy(NoCopy&&) noexcept = default;
	NoCopy& operator=(const NoCopy&) = delete;
	NoCopy& operator=(NoCopy&&) = default;
};

TEST_CASE("zip_inplace") {
	vector<NoCopy<int>> v1;
	v1.emplace_back(1);
	vector<NoCopy<string>> v2;
	v2.emplace_back("I");

	vector<pair<int, string>> result;
	for (auto elem : zip(v1, v2)) {
		result.emplace_back(get<0>(elem).t, get<1>(elem).t);
	}
	vector<pair<int, string>> expected = {{1, "I"}};
	CHECK(result == expected);
}

TEST_CASE("zip_different_lengths") {
	vector<int> v1 = {1, 2, 3, 4, 5};
	const vector<string> v2 = {"I", "II", "III"};

	int count = 0;
	for (auto elem : zip(v1, v2)) {
		count++;
	}
	CHECK(count == 3);

	count = 0;
	for (auto elem : zip(v2, v1)) {
		count++;
	}
	CHECK(count == 3);
}

TEST_CASE("zip_one_empty") {
	vector<int> v1 = {1, 2, 3, 4, 5};
	vector<string> v2;

	int count = 0;
	for (auto elem : zip(v1, v2)) {
		count++;
	}
	CHECK(count == 0);
	for (auto elem : zip(v2, v1)) {
		count++;
	}
	CHECK(count == 0);
}

TEST_CASE("zip_ranges") {
	vector<NoCopy<int>> v1;
	v1.emplace_back(1);

	int count = 0;
	for (auto elem : zip(range(3), range(1, 5))) {
		count++;
	}
	CHECK(count == 3);

	count = 0;
	for (auto elem : zip(v1, range(3))) {
		count++;
	}
	CHECK(count == 1);

	count = 0;
	for (auto elem : zip(range(3), v1)) {
		count++;
	}
	CHECK(count == 1);
}

TEST_CASE("zip_mutable") {
	vector<int> v1 = {1, 2, 3};
	vector<string> v2 = {"I", "II", "III"};

	for (auto elem : zip(v1, v2)) {
		get<0>(elem) = 9;
		get<1>(elem) = "a";
	}
	CHECK(v1 == (vector<int>{9, 9, 9}));
	CHECK(v2 == (vector<string>{"a", "a", "a"}));
}
