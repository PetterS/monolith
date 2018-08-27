#include <list>
#include <map>

#include <catch.hpp>

#include <minimum/core/check.h>
#include <minimum/core/string.h>

using namespace std;
using namespace minimum::core;

void my_failing_function() { minimum_core_assert(1 == 2, "message"); }

TEST_CASE("minimum_core_assert") {
	CHECK_THROWS_AS(minimum_core_assert(1 == 2, "message"), logic_error);
	CHECK_NOTHROW(minimum_core_assert(1 == 1, "message"));

	try {
		my_failing_function();
		FAIL();
	} catch (std::logic_error& err) {
		string message_with_stack_trace = err.what();
		// Will print the stack trace to the test output.
		INFO(message_with_stack_trace);
		CHECK_FALSE(contains(message_with_stack_trace, "verbose_error"));
		CHECK_FALSE(contains(message_with_stack_trace, "get_stack_trace"));
		CHECK(contains(message_with_stack_trace, "my_failing_function"));
		CHECK(contains(message_with_stack_trace, "message"));
	} catch (...) {
		FAIL();
	}
}

TEST_CASE("check") {
	CHECK_THROWS_AS(check(1 == 2, "1 is not 2"), runtime_error);
	CHECK_NOTHROW(check(1 == 1, "Something is very wrong."));
}

TEST_CASE("attest") {
	CHECK_THROWS_AS(attest(1 == 2, "1 is not 2"), logic_error);
	CHECK_NOTHROW(attest(1 == 1, "Something is very wrong."));
}

TEST_CASE("check-message") {
	try {
		check(1 == 2, "1 ", "is ", "not ", 2);
	} catch (std::runtime_error& err) {
		CHECK(err.what() == std::string{"1 is not 2"});
	}
}
