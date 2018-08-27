#include <sstream>
#include <typeinfo>

#include <catch.hpp>

#include <minimum/core/flamegraph.h>

using namespace std;

static void g() {
	FLAMEGRAPH_LOG_FUNCTION;
	{ FLAMEGRAPH_LOG_SCOPE("scope"); }
}

static void f() {
	FLAMEGRAPH_LOG_FUNCTION;
	g();
	{ FLAMEGRAPH_LOG_SCOPE("scope"); }
	{ FLAMEGRAPH_LOG_SCOPE("scope"); }
	g();
}

TEST_CASE("disabled") {
	f();
	CHECK(flamegraph::internal::render_information() == "");
}

TEST_CASE("enabled") {
	flamegraph::enable();
	f();

	auto result = flamegraph::internal::render_information();
	CAPTURE(result);

	CHECK(result.find('f') != string::npos);
	CHECK(result.find('g') != string::npos);
	CHECK(result.find("scope") != string::npos);
}

TEST_CASE("write_pretty") {
	flamegraph::enable();
	f();

	stringstream sout;
	flamegraph::write_pretty(sout);
	auto result = sout.str();
	CAPTURE(result);

	CHECK(result.find('f') != string::npos);
	CHECK(result.find('g') != string::npos);
	CHECK(result.find("scope") != string::npos);
}
