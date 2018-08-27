#include <catch.hpp>

#include <minimum/core/scope_guard.h>

TEST_CASE("scope_exit") {
	int i = 1;
	{ at_scope_exit(i = 3); }
	CHECK(i == 3);

	int j = 1;
	{
		auto guard = ::minimum::core::make_scope_guard([&]() { j = 3; });
		guard.dismiss();
	}
	CHECK(j == 1);

	i = 1;
	j = 1;
	{
		at_scope_exit(i = 3);
		at_scope_exit(j = 3);
		i = 1;
		j = 1;
	}
	CHECK(i == 3);
	CHECK(j == 3);
}
