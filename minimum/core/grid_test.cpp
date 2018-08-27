#include <typeinfo>

#include <catch.hpp>

#include <minimum/core/grid.h>

using namespace std;
using namespace minimum::core;

TEST_CASE("simple_1") {
	auto v = make_grid<double>(10);
	REQUIRE(typeid(v) == typeid(vector<double>));
	REQUIRE(v.size() == 10);
	REQUIRE(v[3] == 0);
}

TEST_CASE("simple_2") {
	auto v = make_grid<double>(10, 20);
	REQUIRE(typeid(v) == typeid(vector<vector<double>>));
	REQUIRE(v.size() == 10);
	REQUIRE(v.back().size() == 20);
	REQUIRE(v[3][3] == 0);
}

TEST_CASE("simple_3") {
	auto v = make_grid<double>(10, 20, 30);
	REQUIRE(typeid(v) == typeid(vector<vector<vector<double>>>));
	REQUIRE(v.size() == 10);
	REQUIRE(v.back().size() == 20);
	REQUIRE(v.back().back().size() == 30);
	REQUIRE(v[3][3][3] == 0);
}

TEST_CASE("lambda_1") {
	int c = 0;
	auto l = [&]() {
		c++;
		return 42.0;
	};
	auto v = make_grid<double>(10, l);
	REQUIRE(c == 10);
	REQUIRE(typeid(v) == typeid(vector<double>));
	REQUIRE(v.size() == 10);
	REQUIRE(v[3] == 42);
}

TEST_CASE("lambda_2") {
	int c = 0;
	auto v = make_grid<double>(10, 20, [&]() {
		c++;
		return 42.0;
	});
	REQUIRE(c == 10 * 20);
	REQUIRE(typeid(v) == typeid(vector<vector<double>>));
	REQUIRE(v.size() == 10);
	REQUIRE(v.back().size() == 20);
	REQUIRE(v[3][3] == 42);
}

TEST_CASE("lambda_3") {
	int c = 0;
	auto v = make_grid<double>(10, 20, 30, [&]() {
		c++;
		return 42.0;
	});
	REQUIRE(c == 10 * 20 * 30);
	REQUIRE(typeid(v) == typeid(vector<vector<vector<double>>>));
	REQUIRE(v.size() == 10);
	REQUIRE(v.back().size() == 20);
	REQUIRE(v.back().back().size() == 30);
	REQUIRE(v[3][3][3] == 42);
}
