
#include <functional>
#include <map>
#include <utility>
#include <vector>

#include <catch.hpp>

#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/ip_to_sat_solver.h>
#include <minimum/linear/pseudoboolean.h>
#include <minimum/linear/pseudoboolean_constraint.h>
#include <minimum/linear/sat_solver.h>
#include <minimum/linear/solver.h>
using minimum::core::contains;
using namespace minimum::linear;

class IPTestFixture {
   public:
	IP ip;
	BooleanVariable x0;
	BooleanVariable x1;
	BooleanVariable x2;
	BooleanVariable x3;
	BooleanVariable x4;
	BooleanVariable x5;

	Sum s0;  // == 0
	Sum s1;  // == 1
	Sum s2;
	Sum s3;
	Sum s4;
	Sum s5;  // == 5

	IPTestFixture() {
		x0 = ip.add_boolean();
		x1 = ip.add_boolean();
		x2 = ip.add_boolean();
		x3 = ip.add_boolean();
		x4 = ip.add_boolean();
		x5 = ip.add_boolean();
		ip.add_constraint(x0 == 1);
		ip.add_constraint(x1 == 1);
		ip.add_constraint(x2 == 1);
		ip.add_constraint(x3 == 1);
		ip.add_constraint(x4 == 1);
		ip.add_constraint(x5 == 1);
		s0 = 0 * x0;
		s1 = 1 * x1;
		s2 = 2 * x2;
		s3 = 3 * x3;
		s4 = 4 * x4;
		s5 = 5 * x5;
		solve(&ip);
	}
};

TEST_CASE_METHOD(IPTestFixture, "addition") {
	PseudoBoolean f = x0 + x1 + x2;
	f -= x1;
	f -= x2;
	CHECK(f.str() == "1*x0 ");
	CHECK((x0 * x1 + x2 * x3).str() == "1*x0*x1 1*x2*x3 ");
	CHECK(((x0 + x1) + x2 * x3).str() == "1*x0 1*x1 1*x2*x3 ");
	CHECK((x2 * x3 + (x0 + x1)).str() == "1*x0 1*x1 1*x2*x3 ");
	CHECK((x1 * x2 + f).str() == "1*x0 1*x1*x2 ");
}

TEST_CASE_METHOD(IPTestFixture, "subtraction") {
	CHECK((x0 * x1 + x0 * x1 - x0 * x2 - x0 * x1 - x0 * x1 + x1).str() == (x1 - x0 * x2).str());

	CHECK((-(x0 * x1) + x1).str() == (x1 - x0 * x1).str());
	PseudoBoolean f = x0 * x1;
	CHECK((-f + x1).str() == (x1 - x0 * x1).str());

	CHECK((x0 * x1 - x2 * x3).str() == "1*x0*x1 -1*x2*x3 ");
	CHECK(((x0 + x1) - x2 * x3).str() == "1*x0 1*x1 -1*x2*x3 ");
	CHECK((x2 * x3 - (x0 + x1)).str() == "-1*x0 -1*x1 1*x2*x3 ");

	CHECK((x0 * x1 - x2).str() == "1*x0*x1 -1*x2 ");
	CHECK((f - f).str() == "");
	CHECK((f - x0 * x1).str() == "");
	CHECK((x0 * x1 - f).str() == "");
}

TEST_CASE("multiplication") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();

	{
		PseudoBoolean f = 2 * x * y * z + x * y;
		CHECK(f.str() == (x * y * z + x * y + x * y * z).str());
		CHECK(f.str() == "1*x0*x1 2*x0*x1*x2 ");
	}

	CHECK((x * x * x).str() == "1*x0 ");

	{
		PseudoBoolean f = x;
		PseudoBoolean g = y;
		CHECK((f + g + 5 * f * g * g * f * f * f).str() == (x + y + 5 * x * y).str());
		CHECK((g + 5 * f * g * g * f * f * f).str() == (y + 5 * x * y).str());
		CHECK((f * g).str() == "1*x0*x1 ");
		CHECK((f * (y * z)).str() == "1*x0*x1*x2 ");
	}

	{
		PseudoBoolean f = x;
		CHECK((2 * (x * y)).str() == "2*x0*x1 ");
		CHECK(((x * y) * 2).str() == "2*x0*x1 ");
		CHECK((2 * f).str() == "2*x0 ");
		CHECK((f * 2).str() == "2*x0 ");
	}

	CHECK(((x + y + 1) * (x + w) * 2).str() == (2 * (x + x * w + y * x + y * w + x + w)).str());
	CHECK(((x + y + 1) * (x + w) * 2).str() == "4*x0 2*x0*x1 2*x0*x3 2*x1*x3 2*x3 ");
	CHECK((x * x * x * x * x * x * x * x * x * x).str() == "1*x0 ");
	CHECK((x * x + x * x + x * x + x * x + x * x).str() == "5*x0 ");
	CHECK((x * y * z * 2).str() == (2 * x * y * z).str());
	CHECK((x * y * 2 * z).str() == (2 * x * y * z).str());
	CHECK((x * 2 * y * z).str() == (2 * x * y * z).str());
	CHECK((x * (2 * y) * z).str() == (2 * x * y * z).str());
	CHECK((((x * 2) * y) * z).str() == (2 * x * y * z).str());
	CHECK(((x * 2) * (y * z)).str() == (2 * x * y * z).str());

	CHECK((x * (y * z)).str() == "1*x0*x1*x2 ");
	CHECK(((x * y) * (z * w)).str() == "1*x0*x1*x2*x3 ");
}

TEST_CASE_METHOD(IPTestFixture, "value") {
	CHECK(PseudoBoolean().value() == 0);
	CHECK(PseudoBoolean(1).value() == 1);
	CHECK(PseudoBoolean(s2).value() == 2);
	CHECK((s1 * s2 * s3 * s4 * s5).value() == 2 * 3 * 4 * 5);
}

TEST_CASE("invalid_type") {
	IP ip;
	auto r = ip.add_variable(IP::Real);
	auto i = ip.add_variable(IP::Integer);
	auto x = ip.add_boolean();
	CHECK_THROWS(x * i);
	CHECK_THROWS(r * x);
	CHECK_THROWS(r * i);
	CHECK_THROWS(x * (x + i));
	CHECK_THROWS(x * x + i);
	CHECK_THROWS(r + x * x);
}

TEST_CASE("identities") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();

	CHECK((x * y + (1 - x) * y + x * (1 - y) + (1 - x) * (1 - y)).str() == "1 ");
	CHECK((x * y * z + x * y * (1 - z) + x * (1 - y) * z + x * (1 - y) * (1 - z) + (1 - x) * y * z
	       + (1 - x) * y * (1 - z) + (1 - x) * (1 - y) * z + (1 - x) * (1 - y) * (1 - z))
	          .str()
	      == "1 ");
}

TEST_CASE("constraints_compilation") {
	IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	auto z = ip.add_boolean();
	auto w = ip.add_boolean();

	ip.add_pseudoboolean_constraint(x * y >= z * w);
	ip.add_pseudoboolean_constraint(x * y >= 3);
	ip.add_pseudoboolean_constraint(3 >= x * y);
	ip.add_pseudoboolean_constraint(x * y >= z + w);
	ip.add_pseudoboolean_constraint(z + w >= x * y);
	ip.add_pseudoboolean_constraint(x * y >= z);
	ip.add_pseudoboolean_constraint(z >= x * y);

	ip.add_pseudoboolean_constraint(x * y <= z * w);
	ip.add_pseudoboolean_constraint(x * y <= 3);
	ip.add_pseudoboolean_constraint(3 <= x * y);
	ip.add_pseudoboolean_constraint(x * y <= z + w);
	ip.add_pseudoboolean_constraint(z + w <= x * y);
	ip.add_pseudoboolean_constraint(x * y <= z);
	ip.add_pseudoboolean_constraint(z <= x * y);

	ip.add_pseudoboolean_constraint(x * y == z * w);
	ip.add_pseudoboolean_constraint(x * y == 3);
	ip.add_pseudoboolean_constraint(3 == x * y);
	ip.add_pseudoboolean_constraint(x * y == z + w);
	ip.add_pseudoboolean_constraint(z + w == x * y);
	ip.add_pseudoboolean_constraint(x * y == z);
	ip.add_pseudoboolean_constraint(z == x * y);

	SUCCEED();
}

TEST_CASE("simple_constraint_0") {
	IP ip;
	ip.add_boolean();
	ip.add_pseudoboolean_constraint(PseudoBoolean{1} == 2);
	try {
		ip.linearize_pseudoboolean_terms();
		FAIL("Should throw.");
	} catch (std::exception& err) {
		CHECK(contains(err.what(), "always false"));
	}
}

TEST_CASE("simple_constraint_1") {
	IP ip;
	auto x = ip.add_boolean();
	ip.add_pseudoboolean_constraint(PseudoBoolean{x} == 1);
	ip.linearize_pseudoboolean_terms();
	solve(&ip);
	CHECK(x.value() == 1);
}

TEST_CASE("simple_constraint_4") {
	IP ip;
	auto x = ip.add_boolean_vector(4);
	ip.add_pseudoboolean_constraint(x[0] * x[1] * x[2] * x[3] == 1);
	ip.linearize_pseudoboolean_terms();
	solve(&ip);
	CHECK(x[0].value() == 1);
	CHECK(x[1].value() == 1);
	CHECK(x[2].value() == 1);
	CHECK(x[3].value() == 1);
	CHECK(ip.get_number_of_variables() == 5);
}

TEST_CASE("simple_objective") {
	IP ip;
	auto x = ip.add_boolean_vector(3);
	ip.add_pseudoboolean_objective(-x[0] * x[1] * x[2]);
	ip.add_pseudoboolean_objective(x[0] * x[1] * x[2]);
	ip.add_pseudoboolean_objective(-x[0] * x[1] * x[2]);
	ip.linearize_pseudoboolean_terms();
	REQUIRE(solve(&ip));
	CHECK(x[0].value() == 1);
	CHECK(x[1].value() == 1);
	CHECK(x[2].value() == 1);
	CHECK(ip.get_number_of_variables() == 4);
}

TEST_CASE("move") {
	IP ip;
	auto x = ip.add_boolean();
	auto constr = x * x * x * x * x >= 0;
	PseudoBooleanConstraint c2 = std::move(constr);
}

TEST_CASE("assignment") {
	PseudoBoolean s1 = Sum(1);
	PseudoBoolean s2 = Sum(2);
	PseudoBoolean s3 = Sum(3);
	s1 = s2;
	s2 = std::move(s3);
	CHECK(s1.value() == 2);
	CHECK(s2.value() == 3);
}

TEST_CASE("golomb") {
	int max_length = 6;  // Maximum length of the ruler.
	int order = 4;       // The order of the created ruler.

	IP ip;
	int num_vars = max_length + 1;
	auto mark = ip.add_boolean_vector(num_vars);
	ip.add_bounds(1, mark[0], 1);
	ip.add_constraint(sum(mark) == order);

	for (int dist = 1; dist < num_vars; ++dist) {
		PseudoBoolean dist_sum;
		for (int i : minimum::core::range(num_vars)) {
			for (int j : minimum::core::range(i + 1, num_vars)) {
				if (j - i == dist) {
					dist_sum += mark[i] * mark[j];
				}
			}
		}
		ip.add_pseudoboolean_constraint(dist_sum <= 1);
	}

	IpToSatSolver solver(bind(minisat_solver, true));
	auto solutions = solver.solutions(&ip);

	int num_solutions = 0;
	while (solutions->get()) {
		num_solutions++;
	}
	CHECK(num_solutions == 2);
}
