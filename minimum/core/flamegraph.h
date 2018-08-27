#pragma once

#include <string>

#include <minimum/core/export.h>

namespace flamegraph {
// Enable flamegraph output for the current thread to the specified filename.
// Different threads should log to different files.
MINIMUM_CORE_API void enable();
// Writes the result for the current thread to the specified filename.
// Stops capturing before writing.
MINIMUM_CORE_API void write_to_file(const std::string& filename);
// Writes the result for the current thread to the specified filename.
// Stops capturing before writing.
MINIMUM_CORE_API void write_pretty(std::ostream& out);
// Place this macro in the beginning of a function to enable logging for the function.
#define FLAMEGRAPH_LOG_FUNCTION flamegraph::internal::Context flamegraph_context(__func__)
// Place this macro in the beginning of a scope to enable logging for that scope.
#define FLAMEGRAPH_LOG_SCOPE(name) flamegraph::internal::Context flamegraph_context(#name)

namespace internal {
MINIMUM_CORE_API bool is_enabled();
MINIMUM_CORE_API void start(const char* name);
MINIMUM_CORE_API void stop();
MINIMUM_CORE_API std::string render_information();

class Context {
   public:
	Context(const char* func) : enabled(is_enabled()) {
		if (enabled) {
			start(func);
		}
	}
	~Context() {
		if (enabled) {
			stop();
		}
	}

   private:
	bool enabled;
};
}  // namespace internal
}  // namespace flamegraph
