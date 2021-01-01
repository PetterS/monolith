// Petter Strandmark
#pragma once

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#include <csignal>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <minimum/core/check.h>
#include <minimum/core/time.h>

namespace minimum {
namespace core {

class Color {
   public:
#ifdef _WIN32
	Color(unsigned short c) : color(c) {}
	unsigned short color;
#else
	Color(const char* s) : str(s) {}
	const char* str;
#endif
};

MINIMUM_CORE_API std::ostream& operator<<(std::ostream& stream, const Color& c);

extern Color MINIMUM_CORE_API NORMAL;
extern const Color MINIMUM_CORE_API WHITE;
extern const Color MINIMUM_CORE_API RED;
extern const Color MINIMUM_CORE_API DKRED;
extern const Color MINIMUM_CORE_API BLUE;
extern const Color MINIMUM_CORE_API DKBLUE;
extern const Color MINIMUM_CORE_API GREEN;
extern const Color MINIMUM_CORE_API DKGREEN;
extern const Color MINIMUM_CORE_API YELLOW;
extern const Color MINIMUM_CORE_API BROWN;

class MINIMUM_CORE_API Timer {
	double start_time = 0;
	double elapsed_time = -1;
	bool active = false;
	void (*old_signal_handler)(int) = SIG_ERR;

   public:
	Timer(const std::string& string);
	~Timer();

	// Stops the timer with a successful message.
	void OK();
	// Prints a message that the operation failed.
	void fail();
	// Ends the current session (if any) and proceeds with
	// the next operation.
	void next(const std::string& string);

	// Obtains the elapsed time after the last success.
	double get_elapsed_time() const;

	// Checks if the operation was interrupted by Ctrl-C.
	// Throws an exception if interrupted.
	void check_for_interruption();
};

template <typename String, typename Function>
void timed_block(const String& string, Function&& function) {
	Timer timer(string);
	function();
	timer.OK();
}
}  // namespace core
}  // namespace minimum
