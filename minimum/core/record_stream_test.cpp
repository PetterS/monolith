#include <sstream>
using namespace std;

#include <catch.hpp>

#include <minimum/core/record_stream.h>
using namespace minimum::core;

TEST_CASE("basic") {
	string empty;
	string small = "Petter";
	string big(500, 'P');

	stringstream stream(ios::binary | ios_base::in | ios_base::out);
	write_record(&stream, small);
	write_record(&stream, small);
	write_record(&stream, big);
	write_record(&stream, empty);
	write_record(&stream, small);
	write_record(&stream, big);

	string tmp;
	read_record(&stream, &tmp);
	CHECK(tmp == small);
	read_record(&stream, &tmp);
	CHECK(tmp == small);
	read_record(&stream, &tmp);
	CHECK(tmp == big);
	read_record(&stream, &tmp);
	CHECK(tmp == empty);
	read_record(&stream, &tmp);
	CHECK(tmp == small);
	read_record(&stream, &tmp);
	CHECK(tmp == big);
}

TEST_CASE("throws_on_read_failure") {
	stringstream stream(ios::binary | ios_base::in | ios_base::out);
	stream << "Petter";
	string tmp;
	CHECK_THROWS(read_record(&stream, &tmp));
}

TEST_CASE("throws_on_write_failure") {
	stringstream stream(ios::binary | ios_base::in);
	CHECK_THROWS(write_record(&stream, "Petter"));
}

class MockProtoBuf {
   public:
	string SerializeAsString() const { return "P"; }
};

TEST_CASE("char_array_correctly_uses_string_function") {
	stringstream stream(ios::binary | ios_base::out);
	write_record(&stream, "Char");
}

TEST_CASE("protobuf") {
	stringstream stream(ios::binary | ios_base::out);
	write_record(&stream, MockProtoBuf{});
}
