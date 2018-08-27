#include <csignal>
#include <utility>

#include <catch.hpp>

#include <minimum/core/death_test.h>

#ifndef WIN32
// These tests work on WIN32, but are too slow (≥1s).
TEST_CASE("abort") { CHECK_DEATH(std::abort()); }

TEST_CASE("terminate") { CHECK_DEATH(std::terminate()); }

TEST_CASE("exit") { CHECK_DEATH(std::exit(0)); }

TEST_CASE("sigsegv") { CHECK_DEATH(std::raise(SIGSEGV)); }
#endif

#ifndef WIN32
// SIGKILL does not exist on WIN32
TEST_CASE("sigkill") { CHECK_DEATH(std::raise(SIGKILL)); }
#endif

TEST_CASE("ok") {
	auto f = []() {};
	CHECK_ALIVE(f());
}

TEST_CASE("exception") { CHECK_ALIVE(throw "Ex"); }

TEST_CASE("exit_code") { CHECK_EXIT_CODE(std::exit(9), 9); }
