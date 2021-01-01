#include <stdexcept>

#include <gflags/gflags.h>
#include <catch.hpp>

#include <minimum/core/main.h>

DEFINE_int32(test_flag_that_is_not_used, 0, "Not used");

TEST_CASE("exception") {
	char* args[] = {const_cast<char*>("prog")};
	CHECK(main_runner([](int, char*[]) -> int { throw "s"; }, 1, args) == 2);
	CHECK(main_runner([](int, char*[]) -> int { throw std::runtime_error("e"); }, 1, args) == 1);
}

TEST_CASE("return_value") {
	char* args[] = {const_cast<char*>("prog")};
	CHECK(main_runner([](int, char*[]) { return 0; }, 1, args) == 0);
	CHECK(main_runner([](int, char*[]) { return 42; }, 1, args) == 42);
}

static int observed_args = -1;
static int observer(int num_args, char*[]) {
	observed_args = num_args;
	return 0;
}

TEST_CASE("gflags_parse") {
	char* args[] = {const_cast<char*>("prog"), const_cast<char*>("--test_flag_that_is_not_used=2")};
	observed_args = -1;
	main_runner(observer, 2, args);
	CHECK(observed_args == 1);
}
