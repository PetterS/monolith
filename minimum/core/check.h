#pragma once

#include <minimum/core/export.h>
#include <minimum/core/string.h>

namespace minimum {
namespace core {

//
// Enables expressions like:
//
//    check(a == 42, a, " is not 42.");
//
// Will throw std::runtime_error if expression is false.
//
template <typename... Args>
void check(bool everything_OK, Args&&... args) {
	if (!everything_OK) {
		throw std::runtime_error(to_string(std::forward<Args>(args)...));
	}
}

void MINIMUM_CORE_API verbose_error_internal(const char* expr,
                                             const char* function_name,
                                             const char* full_file_cstr,
                                             int line,
                                             const char* args);

template <typename... Args>
void verbose_error(const char* expr,
                   const char* function_name,
                   const char* full_file_cstr,
                   int line,
                   Args&&... args) {
	verbose_error_internal(
	    expr, function_name, full_file_cstr, line, to_string(std::forward<Args>(args)...).c_str());
}

// Like check(), but will throw std::logic_error and print a stack trace.
template <typename... Args>
void attest(bool everything_OK, Args&&... args) {
	if (!everything_OK) {
		verbose_error("", "", "", 0, std::forward<Args>(args)...);
	}
}
}  // namespace core
}  // namespace minimum

// Like attest, but will also include the file name, line number, and the expression checked.
#define minimum_core_assert(expr, ...)       \
	(expr) ? ((void)0)                       \
	       : ::minimum::core::verbose_error( \
	             #expr, __func__, __FILE__, __LINE__, ::minimum::core::to_string(__VA_ARGS__))
