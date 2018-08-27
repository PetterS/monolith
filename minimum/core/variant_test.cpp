#include <string_view>
#include <variant>

#include <catch.hpp>

#include <minimum/core/variant.h>

using namespace std;
using namespace minimum::core;

struct S {};

TEST_CASE("overload") {
	variant<string, int> v = 100;
	bool ok = false;
	visit(overload{[&ok](int i) { ok = true; }, [](string_view s) {}}, v);
	CHECK(ok);
}

TEST_CASE("switch_") {
	variant<string, int> v = 100;
	bool ok = false;
	switch_(v,
	        [&ok](int i) {
		        ok = true;
		        CHECK(i == 100);
	        },
	        [](string_view s) { FAIL("Wrong branch called."); });
	CHECK(ok);
}

TEST_CASE("switch_return") {
	variant<string, int> v = "100";
	int result = switch_(v, [](int i) { return 1; }, [](string_view s) { return 2; });
	CHECK(result == 2);
}

constexpr bool test_constexpr() {
	variant<S, int> v = 100;
	int data = 0;
	switch_(v, [&data](int i) { data = i; }, [](S s) {});
	return data == 100;
}

static_assert(test_constexpr(), "constexpr test failed.");
