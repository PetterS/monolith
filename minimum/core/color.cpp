// Petter Strandmark

#include <minimum/core/color.h>

namespace minimum {
namespace core {

std::ostream& operator<<(std::ostream& stream, const Color& c) {
	stream.flush();
#ifdef WIN32
	CONSOLE_SCREEN_BUFFER_INFO buffer_info;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &buffer_info);
	auto background =
	    buffer_info.wAttributes
	    & (BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_BLUE | BACKGROUND_INTENSITY);
	if (NORMAL.color == 0) {
		NORMAL.color =
		    buffer_info.wAttributes
		    & (FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
	}
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), c.color | background);
#else
	stream << "\033[0m" << c.str;
#endif
	stream.flush();

	return stream;
}

#ifdef _WIN32
Color NORMAL = 0;
const Color WHITE = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY;
const Color RED = FOREGROUND_RED | FOREGROUND_INTENSITY;
const Color DKRED = FOREGROUND_RED;
const Color BLUE = FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_INTENSITY;
const Color DKBLUE = FOREGROUND_BLUE | FOREGROUND_GREEN;
const Color GREEN = FOREGROUND_GREEN | FOREGROUND_INTENSITY;
const Color DKGREEN = FOREGROUND_GREEN;
const Color YELLOW = FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY;
const Color BROWN = FOREGROUND_RED | FOREGROUND_GREEN;
#else
Color NORMAL = "";
const Color WHITE = "\033[37;1m";
const Color RED = "\033[31;1m";
const Color DKRED = "\033[31m";
const Color BLUE = "\033[34;1m";
const Color DKBLUE = "\033[34m";
const Color GREEN = "\033[32;1m";
const Color DKGREEN = "\033[32m";
const Color YELLOW = "\033[33;1m";
const Color BROWN = "\033[33m";
#endif

namespace {
volatile std::sig_atomic_t signal_raised = 0;
void interrupt_handler(int signal) { signal_raised = 1; }
}  // namespace

Timer::Timer(const std::string& string) { next(string); }

Timer::~Timer() {
	if (active) {
		fail();
	}
}

void Timer::next(const std::string& string) {
	OK();
	elapsed_time = -1;
	std::cerr << std::left << std::setw(40) << to_string(string, "...") << " [ " << YELLOW << "WAIT"
	          << NORMAL << " ] ";
	start_time = wall_time();

	active = true;
	old_signal_handler = std::signal(SIGINT, interrupt_handler);
	signal_raised = 0;
}

void Timer::OK() {
	if (active) {
		elapsed_time = wall_time() - start_time;
		std::cerr << "\b\b\b\b\b\b\b\b" << GREEN << "  OK  " << NORMAL << "]   ";
		std::cerr << elapsed_time << " s." << std::endl;
		active = false;
		if (old_signal_handler != SIG_ERR) {
			std::signal(SIGINT, old_signal_handler);
		}
	}
}

double Timer::get_elapsed_time() const {
	minimum_core_assert(elapsed_time >= 0);
	return elapsed_time;
}

void Timer::fail() {
	if (active) {
		std::cerr << "\b\b\b\b\b\b\b\b" << RED << "FAILED" << NORMAL << "]" << std::endl;
		active = false;
		if (old_signal_handler != SIG_ERR) {
			std::signal(SIGINT, old_signal_handler);
		}
	}
}

void Timer::check_for_interruption() {
	if (signal_raised == 1 && active) {
		std::cerr << "\b\b\b\b\b\b\b\b" << RED << "Ctrl+C" << NORMAL << "]" << std::endl;
		active = false;
		if (old_signal_handler != SIG_ERR) {
			std::signal(SIGINT, old_signal_handler);
		}
		throw std::runtime_error("Operation interrupted.");
	}
}
}  // namespace core
}  // namespace minimum
