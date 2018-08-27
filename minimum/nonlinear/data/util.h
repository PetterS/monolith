#pragma once

#include <string>

#include <minimum/core/check.h>

namespace minimum {
namespace nonlinear {
namespace data {
namespace internal {
namespace {
std::string get_directory(std::string this_file_name) {
	auto sep = this_file_name.find_last_of("/\\");
	minimum_core_assert(sep != std::string::npos);
	return this_file_name.substr(0, sep);
}
}  // namespace
}  // namespace internal
namespace {
std::string get_directory() { return internal::get_directory(__FILE__); }
}  // namespace
}  // namespace data
}  // namespace nonlinear
}  // namespace minimum
