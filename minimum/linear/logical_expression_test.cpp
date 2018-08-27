// Petter Strandmark 2013.

#include <catch.hpp>

#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>

using namespace std;
using namespace minimum::linear;

class IPTestFixture {
   public:
	IP ip;
	BooleanVariable x0;  // == 0
	BooleanVariable x1;  // == 1
	BooleanVariable y0;  // == 0
	BooleanVariable y1;  // == 1

	IPTestFixture()
	    : x0(ip.add_boolean()), x1(ip.add_boolean()), y0(ip.add_boolean()), y1(ip.add_boolean()) {
		ip.add_bounds(0, x0, 0);
		ip.add_bounds(1, x1, 1);
		ip.add_bounds(0, y0, 0);
		ip.add_bounds(1, y1, 1);
		ip.add_constraint(x0 || x1 || y0 || y1);
		solve(&ip);
	}
};

TEST_CASE_METHOD(IPTestFixture, "test_variables") {
	CHECK(ip.get_solution(x0) == false);
	CHECK(ip.get_solution(x1) == true);
	CHECK(ip.get_solution(y0) == false);
	CHECK(ip.get_solution(y1) == true);
}

TEST_CASE_METHOD(IPTestFixture, "test_not") {
	CHECK(ip.get_solution(!x0) == true);
	CHECK(ip.get_solution(!x1) == false);
}

TEST_CASE_METHOD(IPTestFixture, "test_or") {
	CHECK(ip.get_solution(x1 || y1) == true);
	CHECK(ip.get_solution(x1 || x1) == true);
	CHECK(ip.get_solution(x1 || x1 || y1) == true);
	CHECK(ip.get_solution(x1 || y1 || x1) == true);
	CHECK(ip.get_solution(y1 || x1 || x1) == true);
	CHECK(ip.get_solution(x0 || x1 || y0) == true);
	CHECK(ip.get_solution(x1 || x0 || y0) == true);
	CHECK(ip.get_solution(y0 || x0 || y1) == true);

	CHECK(ip.get_solution(x0 || y0) == false);
	CHECK(ip.get_solution(y0 || x0 || y0) == false);
	CHECK(ip.get_solution(x0 || x0 || y0) == false);
}

TEST_CASE_METHOD(IPTestFixture, "test_or_overloads") {
	auto s0 = x0 || y0;
	auto s1 = x1 || y1;
	CHECK(ip.get_solution(s0 || s1) == true);
	CHECK(ip.get_solution(s0 || x1) == true);
	CHECK(ip.get_solution(x0 || s1) == true);
}

TEST_CASE_METHOD(IPTestFixture, "implication") {
	auto s0 = x0 || y0;
	CHECK(ip.get_solution(implication(x1, s0)) == false);
	CHECK(ip.get_solution(implication(x1, y1 || x0)) == true);
}
