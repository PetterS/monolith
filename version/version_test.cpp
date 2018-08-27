#include <catch.hpp>

#include <version/version.h>

TEST_CASE("Version") {
	CAPTURE(version::revision);
	CHECK(!version::revision.empty());
}
