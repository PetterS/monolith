#include <list>
#include <map>

#include <catch.hpp>

#include <minimum/core/color.h>

TEST_CASE("timer") {
	minimum::core::Timer t("Testing timer");
	CHECK_THROWS(t.get_elapsed_time());
	t.OK();
	CHECK(t.get_elapsed_time() >= 0);
}

TEST_CASE("timer_explicit_fail") {
	minimum::core::Timer t("Testing timer");
	t.fail();
	CHECK_THROWS(t.get_elapsed_time());
}

TEST_CASE("timer_auto_fail") {
	CHECK_NOTHROW([]() { minimum::core::Timer t("Testing timer"); }());
}

TEST_CASE("timed_block") {
	minimum::core::timed_block("Testing block.", []() {});
}

static void empty_handler(int) {}

TEST_CASE("signal") {
	std::signal(SIGINT, empty_handler);
	REQUIRE((std::signal(SIGINT, empty_handler) == empty_handler));
	{
		minimum::core::Timer t("Testing timer");
		CHECK_NOTHROW(t.check_for_interruption());
		std::raise(SIGINT);
		CHECK_THROWS(t.check_for_interruption());
		t.OK();
	}
	REQUIRE((std::signal(SIGINT, SIG_DFL) == empty_handler));
}
