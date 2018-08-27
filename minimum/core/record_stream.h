#pragma once
#include <iostream>
#include <string>
#include <type_traits>

#include <minimum/core/check.h>
#include <minimum/core/export.h>

namespace minimum {
namespace core {
MINIMUM_CORE_API void write_record(std::ostream* out, const std::string& data);
MINIMUM_CORE_API void read_record(std::istream* in, std::string* data);

// Reads a protobuf message from the record stream.
template <typename Proto>
Proto read_record(std::istream* in) {
	Proto proto;
	std::string data;
	read_record(in, &data);
	check(proto.ParseFromString(data), "Failed to parse proto.");
	return proto;
}

// Writes a protobuf message to the record stream.
template <typename Proto>
std::void_t<decltype(&Proto::SerializeAsString)> write_record(std::ostream* out,
                                                              const Proto& proto) {
	write_record(out, proto.SerializeAsString());
}
}  // namespace core
}  // namespace minimum
