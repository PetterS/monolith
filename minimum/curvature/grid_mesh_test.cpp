// Petter Strandmark 2013.
#include <algorithm>
#include <cmath>

#include <catch.hpp>

#include <minimum/curvature/grid_mesh.h>

using namespace minimum::curvature;

void add_edges_test(float distance, int expected_connectivity) {
	const int n = 10;
	GridMesh mesh(n, n, distance, false);
	CHECK(mesh.get_connectivity() == expected_connectivity);
}

TEST_CASE("GridMesh/add_edges", "") {
	add_edges_test(std::sqrt(2.0f * 2.0f + 1.0f * 1.0f), 16);
	add_edges_test(std::sqrt(2.0f), 8);
	add_edges_test(1.0f, 4);
}
