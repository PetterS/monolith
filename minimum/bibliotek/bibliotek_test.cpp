
#include <catch.hpp>

#include <minimum/bibliotek/bibliotek.h>
#include <minimum/linear/ip.h>

TEST_CASE("ip") {
	minimum::linear::IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x + y >= 1);

	auto ip_str = ip.get().SerializeAsString();
	auto ip_ptr = minimum_linear_ip(ip_str.data(), ip_str.size());
	CHECK(ip_ptr != nullptr);
	minimum_linear_ip_destroy(ip_ptr);
}

TEST_CASE("solver") {
	auto solver_ptr = minimum_linear_solver_default();
	CHECK(solver_ptr != nullptr);
	minimum_linear_solver_destroy(solver_ptr);
}

TEST_CASE("solutions") {
	minimum::linear::IP ip;
	auto x = ip.add_boolean();
	auto y = ip.add_boolean();
	ip.add_constraint(x + y >= 1);

	auto ip_str = ip.get().SerializeAsString();
	auto ip_ptr = minimum_linear_ip(ip_str.data(), ip_str.size());
	auto solver_ptr = minimum_linear_solver_default();
	auto solutions_ptr = minimum_linear_solutions(solver_ptr, ip_ptr);

	std::vector<double> solution(2, 0.0);
	auto result =
	    minimum_linear_solutions_get(solutions_ptr, ip_ptr, solution.data(), solution.size());
	CHECK(result == 0);

	minimum_linear_solutions_destroy(solutions_ptr);
	minimum_linear_solver_destroy(solver_ptr);
	minimum_linear_ip_destroy(ip_ptr);
}
