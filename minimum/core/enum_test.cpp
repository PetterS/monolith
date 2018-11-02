#include <iostream>

#include <gflags/gflags.h>
#include <catch.hpp>

#include <minimum/core/enum.h>

MAKE_ENUM(MyEnum, FOO, BAR, BAZ, BIZ);

TEST_CASE("to_string") {
	CHECK(to_string(MyEnum::FOO) == "FOO");
	CHECK(to_string(MyEnum::BAZ) == "BAZ");
}

TEST_CASE("from_string") {
	CHECK(MyEnum_from_string("BAR") == MyEnum::BAR);
	CHECK(MyEnum_from_string("BIZ") == MyEnum::BIZ);
	CHECK_THROWS(MyEnum_from_string("Petter"));
}

TEST_CASE("from_int") {
	CHECK(MyEnum_from_int(1) == MyEnum::BAR);
	CHECK_THROWS(MyEnum_from_int(12));
	MyEnum_enum_names;
}

DEFINE_enum_flag(MyEnum, my_flag, MyEnum::BAR, "My enum flag 2.");

TEST_CASE("flag") {
	CHECK(gflags::ReadFlagsFromString("--my_flag=BAZ", "test", false));
	CHECK_FALSE(gflags::ReadFlagsFromString("--my_flag=Petter", "test", false));
	CHECK(gflags::ReadFlagsFromString("--my_flag=BIZ", "test", false));
	CHECK(FLAGS_my_flag == "BIZ");
	CHECK(MyEnum_from_string(FLAGS_my_flag) == MyEnum::BIZ);
}

MAKE_ENUM(MyEnum1, A1);
MAKE_ENUM(MyEnum2, A1, A2);
MAKE_ENUM(MyEnum9, A1, A2, A3, A4, A5, A6, A7, A8, A9);

// To disable warnings.
void use_all_functions() {
#define USE(Name)           \
	to_string(Name{});      \
	Name##_from_string(""); \
	Name##_from_int(0);     \
	Name##_flag_validator(nullptr, "");

	USE(MyEnum1);
	USE(MyEnum2);
	USE(MyEnum9);
}
