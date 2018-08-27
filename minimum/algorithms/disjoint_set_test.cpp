#include <vector>
using namespace std;

#include <catch.hpp>

#include <minimum/algorithms/disjoint_set.h>
#include <minimum/core/range.h>
using namespace minimum::algorithms;
using namespace minimum::core;

TEST_CASE("self_as_parent") {
	DisjointSetElement a;
	CHECK(a.find() == &a);
}

TEST_CASE("merge2") {
	DisjointSetElement a, b;
	CHECK(a.join(&b));
	CHECK_FALSE(a.join(&b));
	CHECK(a.find() == b.find());
}

TEST_CASE("10") {
	int n = 10;
	vector<DisjointSetElement> elem(n);
	for (auto i : range(n)) {
		CHECK(elem[i].find() == &elem[i]);
	}

	CHECK(elem[0].join(&elem[1]));
	CHECK(elem[1].join(&elem[2]));
	CHECK(elem[2].join(&elem[3]));

	CHECK(elem[4].join(&elem[5]));
	CHECK_FALSE(elem[5].join(&elem[4]));
	CHECK(elem[6].join(&elem[4]));

	CHECK(elem[5].join(&elem[1]));

	CHECK(elem[7].join(&elem[9]));
	CHECK_FALSE(elem[8].join(&elem[8]));
	CHECK(elem[9].join(&elem[8]));

	CHECK(elem[0].join(&elem[9]));

	auto parent = elem[0].find();
	for (auto i : range(n)) {
		CHECK(elem[i].find() == parent);
	}
}
