#include <catch.hpp>

#include <minimum/core/read_matrix.h>
#include <minimum/core/time.h>
#include <minimum/core/write_matrix.h>
#include <iostream>

using namespace std;
using minimum::core::read_matrix;
using minimum::core::wall_time;
using minimum::core::write_matrix;
using std::vector;

TEST_CASE("small") {
	stringstream ss;
	vector<vector<int>> mat{{1, 2}, {3, 4}};
	write_matrix(ss, mat);
	CHECK(read_matrix<int>(ss) == mat);
}

TEST_CASE("large") {
	int n = 1021;
	int m = 768;
	stringstream ss;
	vector<vector<int>> mat;
	for (int i = 0; i < m; ++i) {
		mat.emplace_back();
		for (int j = 0; j < n; ++j) {
			mat.back().emplace_back(i * j);
		}
	}
	double start_time = wall_time();
	write_matrix(ss, mat);
	CAPTURE(wall_time() - start_time);
	SUCCEED();

	start_time = wall_time();
	auto res = read_matrix<int>(ss);
	CAPTURE(wall_time() - start_time);
	CHECK((res == mat));
}
