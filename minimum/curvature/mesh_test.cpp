// Petter Strandmark 2013.
#include <algorithm>
#include <cmath>

#include <catch.hpp>

#include <minimum/curvature/mesh.h>

using namespace minimum::curvature;

#define EXPECT_EQ(a, b) CHECK((a) == (b))
#define EXPECT_FLOAT_EQ(a, b) CHECK(Approx(float(a)) == float(b))
#define EXPECT_LT(a, b) CHECK((a) < (b))
#define EXPECT_NE(a, b) CHECK((a) != (b))

#define ASSERT_EQ(a, b) REQUIRE((a) == (b))
#define ASSERT_GE(a, b) REQUIRE((a) >= (b))
#define ASSERT_NE(a, b) REQUIRE((a) != (b))

TEST_CASE("Mesh/find_point", "") {
	Mesh mesh;
	const int n = 10;

	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			for (int z = 0; z < n; ++z) {
				mesh.add_point(x, y, z);
			}
		}
	}

	EXPECT_EQ(mesh.number_of_points(), n * n * n);

	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			for (int z = 0; z < n; ++z) {
				// Find the exact point.
				int i = mesh.find_point(x, y, z);
				ASSERT_GE(i, 0);
				EXPECT_EQ(mesh.get_point(i).x, x);
				EXPECT_EQ(mesh.get_point(i).y, y);
				EXPECT_EQ(mesh.get_point(i).z, z);

				// Search for nonexistant points.
				i = mesh.find_point(x + 0.5, y, z);
				EXPECT_LT(i, 0);
				i = mesh.find_point(x, y + 0.5, z);
				EXPECT_LT(i, 0);
				i = mesh.find_point(x, y, z + 0.5);
				EXPECT_LT(i, 0);

				// Search for slightly perturbed versions
				// of the point.
				i = mesh.find_point(x + 1e-5, y, z);
				ASSERT_GE(i, 0);
				EXPECT_EQ(mesh.get_point(i).x, x);
				EXPECT_EQ(mesh.get_point(i).y, y);
				EXPECT_EQ(mesh.get_point(i).z, z);

				i = mesh.find_point(x, y + 1e-5, z);
				ASSERT_GE(i, 0);
				EXPECT_EQ(mesh.get_point(i).x, x);
				EXPECT_EQ(mesh.get_point(i).y, y);
				EXPECT_EQ(mesh.get_point(i).z, z);

				i = mesh.find_point(x + 1e-5, y, z + 1e-5);
				ASSERT_GE(i, 0);
				EXPECT_EQ(mesh.get_point(i).x, x);
				EXPECT_EQ(mesh.get_point(i).y, y);
				EXPECT_EQ(mesh.get_point(i).z, z);
			}
		}
	}
}

TEST_CASE("Mesh/add_edge", "") {
	Mesh mesh;
	mesh.add_point(0, 0, 0);
	mesh.add_point(1, 0, 0);
	mesh.add_point(0, 1, 0);

	mesh.add_edge(0, 0, 0, 1, 0, 0);
	mesh.add_edge(0, 0, 0, 0, 1, 0);
	mesh.add_edge(1, 0, 0, 0, 1, 0);
	mesh.add_edge(1, 0, 0, 0, 1, 0);
	mesh.add_edge(0, 1, 0, 1, 0, 0);
	mesh.add_edge(1, 0, 0, 0, 1, 0.00001f);

	// Adding an edge from p to p is an error.
	CHECK_THROWS_AS(mesh.add_edge(1, 0, 0, 1, 0, 0), std::runtime_error);

	EXPECT_EQ(mesh.number_of_edges(), 2 * 3);
}

TEST_CASE("Mesh/edge_length", "") {
	Mesh mesh;
	mesh.add_point(0, 0, 0);
	mesh.add_point(1, 1, 0);
	mesh.add_edge(0, 0, 0, 1, 1, 0);
	EXPECT_FLOAT_EQ(mesh.edge_length(0), std::sqrt(2.0f));
}

TEST_CASE("Mesh/create_edge_pairs3", "") {
	Mesh mesh;
	mesh.add_point(0, 0, 0);
	mesh.add_point(1, 0, 0);
	mesh.add_point(0, 1, 0);

	mesh.add_edge(0, 0, 0, 1, 0, 0);
	mesh.add_edge(0, 0, 0, 0, 1, 0);
	mesh.add_edge(1, 0, 0, 0, 1, 0);

	mesh.finish();

	ASSERT_EQ(mesh.number_of_edges(), 2 * 3);

	std::vector<int> adjacent;
	mesh.get_adjacent_edges(0, &adjacent);
	EXPECT_EQ(adjacent.size(), 1);

	mesh.get_adjacent_edges(1, &adjacent);
	EXPECT_EQ(adjacent.size(), 1);

	mesh.get_adjacent_edges(2, &adjacent);
	EXPECT_EQ(adjacent.size(), 1);
}

class MeshLine : public Mesh {
   public:
	virtual void SetUp() {
		add_point(0, 0, 0);
		add_point(1, 0, 0);
		add_point(2, 0, 0);
		add_point(3, 0, 0);
		add_point(4, 0, 0);

		add_edge(0, 0, 0, 1, 0, 0);
		add_edge(1, 0, 0, 2, 0, 0);
		add_edge(2, 0, 0, 3, 0, 0);
		add_edge(3, 0, 0, 4, 0, 0);

		finish(true);
	}
};

TEST_CASE("MeshLine/number_of_edges", "") {
	MeshLine mesh;
	mesh.SetUp();

	EXPECT_EQ(mesh.number_of_edges(), 2 * 4);
	EXPECT_EQ(mesh.get_connectivity(), 2);
}

TEST_CASE("MeshLine/number_of_edge_neighbors", "") {
	MeshLine mesh;
	mesh.SetUp();

	std::vector<int> adjacent;
	int e = mesh.find_edge(0, 1);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 1);

	e = mesh.find_edge(1, 0);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 0);

	e = mesh.find_edge(1, 2);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 1);

	e = mesh.find_edge(2, 1);
	EXPECT_EQ(adjacent.size(), 1);
}

TEST_CASE("MeshLine/number_of_edge_pairs", "") {
	MeshLine mesh;
	mesh.SetUp();

	// for (int ep = 0; ep < mesh.number_of_edge_pairs(); ++ep) {
	//	int p1 = mesh.get_edge_pair(ep).p1;
	//	int p2 = mesh.get_edge_pair(ep).p2;
	//	int p3 = mesh.get_edge_pair(ep).p3;
	//	std::cerr << "(" << p1 << ", " << p2 << ", " << p3 << ")\n";
	//}
	EXPECT_EQ(mesh.number_of_edge_pairs(), 6);
}

TEST_CASE("MeshLine/edgepath_to_points", "") {
	MeshLine mesh;
	mesh.SetUp();

	int p0 = mesh.find_point(0, 0, 0);
	int p1 = mesh.find_point(1, 0, 0);
	int p2 = mesh.find_point(2, 0, 0);
	int p3 = mesh.find_point(3, 0, 0);
	int p4 = mesh.find_point(4, 0, 0);

	std::vector<int> edges(4);
	edges[0] = mesh.find_edge(p0, p1);
	edges[1] = mesh.find_edge(p1, p2);
	edges[2] = mesh.find_edge(p2, p3);
	edges[3] = mesh.find_edge(p3, p4);

	auto points = mesh.edgepath_to_points(edges);
	ASSERT_EQ(points.size(), 5);
	EXPECT_FLOAT_EQ(points[0].x, 0);
	EXPECT_FLOAT_EQ(points[1].x, 1);
	EXPECT_FLOAT_EQ(points[2].x, 2);
	EXPECT_FLOAT_EQ(points[3].x, 3);
	EXPECT_FLOAT_EQ(points[4].x, 4);
}

TEST_CASE("MeshLine/pairpath_to_points", "") {
	MeshLine mesh;
	mesh.SetUp();

	int p0 = mesh.find_point(0, 0, 0);
	int p1 = mesh.find_point(1, 0, 0);
	int p2 = mesh.find_point(2, 0, 0);
	int p3 = mesh.find_point(3, 0, 0);
	int p4 = mesh.find_point(4, 0, 0);

	std::vector<int> pairs(3);
	pairs[0] = mesh.find_edge_pair(p0, p1, p2);
	pairs[1] = mesh.find_edge_pair(p1, p2, p3);
	pairs[2] = mesh.find_edge_pair(p2, p3, p4);

	auto points = mesh.pairpath_to_points(pairs);
	ASSERT_EQ(points.size(), 5);
	EXPECT_FLOAT_EQ(points[0].x, 0);
	EXPECT_FLOAT_EQ(points[1].x, 1);
	EXPECT_FLOAT_EQ(points[2].x, 2);
	EXPECT_FLOAT_EQ(points[3].x, 3);
	EXPECT_FLOAT_EQ(points[4].x, 4);
}

class Mesh4 : public Mesh {
   public:
	Mesh4() {
		//
		//  0---1
		//  |\ /|
		//  | X |
		//  |/ \|
		//  2---3
		//
		add_point(0, 0, 0);
		add_point(1, 0, 0);
		add_point(0, 1, 0);
		add_point(1, 1, 0);

		add_edge(0, 0, 0, 1, 0, 0);
		add_edge(0, 0, 0, 0, 1, 0);
		add_edge(1, 0, 0, 1, 1, 0);
		add_edge(0, 1, 0, 1, 1, 0);
		add_edge(0, 0, 0, 1, 1, 0);
		add_edge(0, 1, 0, 1, 0, 0);

		finish(true);
	}
};

TEST_CASE("Mesh4/finish", "") {
	Mesh4 mesh;

	EXPECT_EQ(mesh.number_of_points(), 4);
	EXPECT_EQ(mesh.number_of_edges(), 2 * 6);
}

TEST_CASE("Mesh4/number_of_edge_neighbors", "") {
	Mesh4 mesh;

	std::vector<int> adjacent;
	int e;

	auto check_contains = [](const std::vector<int>& v, int value) {
		auto itr = std::find(v.begin(), v.end(), value);
		REQUIRE(itr != v.end());
		EXPECT_EQ(*itr, value);
	};

	e = mesh.find_edge(0, 1);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 2);
	check_contains(adjacent, mesh.find_edge(1, 3));
	check_contains(adjacent, mesh.find_edge(1, 2));

	e = mesh.find_edge(1, 0);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 2);
	check_contains(adjacent, mesh.find_edge(0, 3));
	check_contains(adjacent, mesh.find_edge(0, 2));

	e = mesh.find_edge(0, 3);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 2);
	check_contains(adjacent, mesh.find_edge(3, 1));
	check_contains(adjacent, mesh.find_edge(3, 2));

	e = mesh.find_edge(2, 1);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 2);
	check_contains(adjacent, mesh.find_edge(1, 0));
	check_contains(adjacent, mesh.find_edge(1, 3));

	e = mesh.find_edge(1, 3);
	mesh.get_adjacent_edges(e, &adjacent);
	EXPECT_EQ(adjacent.size(), 2);
	check_contains(adjacent, mesh.find_edge(3, 0));
	check_contains(adjacent, mesh.find_edge(3, 2));
}

TEST_CASE("Mesh4/edge_adjacency", "") {
	Mesh4 mesh;

	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		std::vector<int> adjacency;
		mesh.get_adjacent_edges(e, &adjacency);

		for (auto itr = adjacency.begin(); itr != adjacency.end(); ++itr) {
			int p1 = mesh.get_edge(e).first;
			int p2 = mesh.get_edge(e).second;
			int p3 = mesh.get_edge(*itr).first;
			int p4 = mesh.get_edge(*itr).second;
			ASSERT_NE(p1, p2);
			ASSERT_NE(p1, p3);
			ASSERT_NE(p1, p4);
			ASSERT_EQ(p2, p3);
			ASSERT_NE(p3, p4);
		}
	}
}

TEST_CASE("Mesh4/number_of_pair_neighbors", "") {
	Mesh4 mesh;

	std::vector<int> adjacent;
	int e;

	auto check_contains = [](const std::vector<int>& v, int value) {
		auto itr = std::find(v.begin(), v.end(), value);
		REQUIRE(itr != v.end());
		EXPECT_EQ(*itr, value);
	};

	e = mesh.find_edge_pair(0, 1, 3);
	mesh.get_adjacent_pairs(e, &adjacent);
	// for (int ep : adjacent) {
	//	std::cerr << "("
	//	          << std::get<0>(mesh.get_edge_pair(ep)) << ','
	//	          << std::get<1>(mesh.get_edge_pair(ep)) << ','
	//	          << std::get<2>(mesh.get_edge_pair(ep)) << ")\n";
	//}
	// Only one neighboring pair, because (3, 2, 0) overlaps
	// with (0, 1, 3).
	EXPECT_EQ(adjacent.size(), 1);
	check_contains(adjacent, mesh.find_edge_pair(1, 3, 2));
}

TEST_CASE("Mesh4/pair_adjacency", "") {
	Mesh4 mesh;

	using std::tie;

	std::vector<int> adjacency;
	adjacency.reserve(100);

	for (int ep = 0; ep < mesh.number_of_edge_pairs(); ++ep) {
		mesh.get_adjacent_pairs(ep, &adjacency);

		for (auto itr = adjacency.begin(); itr != adjacency.end(); ++itr) {
			int p1, p2, p3, p4, p5, p6;
			tie(p1, p2, p3) = mesh.get_edge_pair(ep);
			tie(p4, p5, p6) = mesh.get_edge_pair(*itr);

			EXPECT_NE(p1, p2);
			EXPECT_NE(p1, p3);
			EXPECT_NE(p1, p4);
			EXPECT_NE(p1, p5);
			EXPECT_NE(p1, p6);

			EXPECT_NE(p6, p1);
			EXPECT_NE(p6, p2);
			EXPECT_NE(p6, p3);
			EXPECT_NE(p6, p4);
			EXPECT_NE(p6, p5);

			EXPECT_EQ(p2, p4);
			EXPECT_EQ(p3, p5);
		}
	}
}

class MeshComplete : public Mesh {
   public:
	int n;

	MeshComplete() {
		n = 8;
		for (int i = 0; i < n; ++i) {
			add_point(i, 0, 0);
		}

		for (int i = 0; i < n; ++i) {
			for (int j = i + 1; j < n; ++j) {
				add_edge(i, 0, 0, j, 0, 0);
			}
		}
		finish(true);
	}
};

TEST_CASE("MeshComplete/finish", "") {
	MeshComplete mesh;
	int n = mesh.n;

	EXPECT_EQ(mesh.number_of_points(), n);
	EXPECT_EQ(mesh.number_of_edges(), n * (n - 1));
	EXPECT_EQ(mesh.get_connectivity(), n - 1);
}

TEST_CASE("MeshComplete/edge_adjacency", "") {
	MeshComplete mesh;
	int n = mesh.n;

	for (int e = 0; e < mesh.number_of_edges(); ++e) {
		std::vector<int> adjacency;
		mesh.get_adjacent_edges(e, &adjacency);

		for (auto itr = adjacency.begin(); itr != adjacency.end(); ++itr) {
			int p1 = mesh.get_edge(e).first;
			int p2 = mesh.get_edge(e).second;
			int p3 = mesh.get_edge(*itr).first;
			int p4 = mesh.get_edge(*itr).second;
			ASSERT_NE(p1, p2);
			ASSERT_NE(p1, p3);
			ASSERT_NE(p1, p4);
			ASSERT_EQ(p2, p3);
			ASSERT_NE(p3, p4);
		}
	}
}

TEST_CASE("MeshComplete/number_of_pair_neighbors", "") {
	MeshComplete mesh;
	int n = mesh.n;

	std::vector<int> adjacent;

	auto check_contains = [](const std::vector<int>& v, int value) {
		auto itr = std::find(v.begin(), v.end(), value);
		REQUIRE(itr != v.end());
		EXPECT_EQ(*itr, value);
	};

	int e = mesh.find_edge_pair(0, 1, 2);
	mesh.get_adjacent_pairs(e, &adjacent);
	EXPECT_EQ(adjacent.size(), n - 3);
	for (int i = 3; i < n; ++i) {
		check_contains(adjacent, mesh.find_edge_pair(1, 2, i));
	}

	for (int p1 = 0; p1 < n; ++p1) {
		for (int p2 = 0; p2 < n; ++p2) {
			for (int p3 = 0; p3 < n; ++p3) {
				if (p1 != p2 && p1 != p3 && p2 != p3) {
					int e = mesh.find_edge_pair(p1, p2, p3);
					mesh.get_adjacent_pairs(e, &adjacent);
					EXPECT_EQ(adjacent.size(), n - 3);
					for (int p4 = 0; p4 < n; ++p4) {
						if (p4 != p1 && p4 != p2 && p4 != p3) {
							check_contains(adjacent, mesh.find_edge_pair(p2, p3, p4));
						}
					}
				}
			}
		}
	}
}

TEST_CASE("MeshComplete/pair_adjacency", "") {
	MeshComplete mesh;
	int n = mesh.n;

	using std::tie;

	std::vector<int> adjacency;
	adjacency.reserve(100);

	for (int ep = 0; ep < mesh.number_of_edge_pairs(); ++ep) {
		mesh.get_adjacent_pairs(ep, &adjacency);

		for (auto itr = adjacency.begin(); itr != adjacency.end(); ++itr) {
			int p1, p2, p3, p4, p5, p6;
			tie(p1, p2, p3) = mesh.get_edge_pair(ep);
			tie(p4, p5, p6) = mesh.get_edge_pair(*itr);

			EXPECT_NE(p1, p2);
			EXPECT_NE(p1, p3);
			EXPECT_NE(p1, p4);
			EXPECT_NE(p1, p5);
			EXPECT_NE(p1, p6);

			EXPECT_NE(p6, p1);
			EXPECT_NE(p6, p2);
			EXPECT_NE(p6, p3);
			EXPECT_NE(p6, p4);
			EXPECT_NE(p6, p5);

			EXPECT_EQ(p2, p4);
			EXPECT_EQ(p3, p5);
		}
	}
}

void add_edges_test(float distance, int expected_connectivity) {
	int n = 10;
	Mesh mesh;
	for (int x = 0; x < n; ++x) {
		for (int y = 0; y < n; ++y) {
			mesh.add_point(float(x), float(y), 0.0f);
		}
	}
	mesh.add_edges(distance + 0.1f);
	mesh.finish();

	EXPECT_EQ(mesh.get_connectivity(), expected_connectivity);
}

TEST_CASE("Mesh/add_edges", "") {
	add_edges_test(std::sqrt(2.0f * 2.0f + 1.0f * 1.0f), 16);
	add_edges_test(std::sqrt(2.0f), 8);
	add_edges_test(1.0f, 4);
}
