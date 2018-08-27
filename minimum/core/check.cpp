#include <cstdlib>
#include <ctime>
#include <mutex>

#include <minimum/core/check.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
// clang-format off
#include <windows.h>
#include <Dbghelp.h>
// clang-format on
#pragma comment(lib, "Dbghelp.lib")
#elif !defined(__CYGWIN__) && !defined(EMSCRIPTEN)
#include <cxxabi.h>
#include <execinfo.h>
#endif

namespace minimum {
namespace core {

std::mutex minimum_core_stack_trace_mutex;

std::string get_stack_trace(void) {
	using namespace std;

	// Absolutely required for the Windows version.
	lock_guard<mutex> lock(minimum_core_stack_trace_mutex);

	const unsigned int max_stack_size = 10000;
	void* stack[max_stack_size];

	string sout = "\n\nStack trace:\n";

#ifdef _WIN32
	HANDLE process;
	process = GetCurrentProcess();
	SymInitialize(process, NULL, TRUE);

	unsigned short frames = CaptureStackBackTrace(0, max_stack_size, stack, NULL);

	typedef IMAGEHLP_SYMBOL64 SymbolInfo;
	SymbolInfo* symbol = (SymbolInfo*)malloc(sizeof(SymbolInfo) + 256);
	symbol->MaxNameLength = 255;
	symbol->SizeOfStruct = sizeof(SymbolInfo);

	for (unsigned int i = 0; i < frames; i++) {
		string symbol_name = "???";
		ptrdiff_t address = ptrdiff_t(stack[i]);
		ptrdiff_t offset = 0;
		if (SymGetSymFromAddr64(process, DWORD64(stack[i]), nullptr, symbol)) {
			address = ptrdiff_t(symbol->Address);
			symbol_name = symbol->Name;
			offset = ptrdiff_t(stack[i]) - address;
		}

		if (symbol_name.find("::get_stack_trace") != symbol_name.npos) {
			continue;
		}
		if (symbol_name.find("::verbose_error") != symbol_name.npos) {
			continue;
		}

		sout += to_string(frames - i - 1, ": ", symbol_name, " + ", offset, "\n");
	}
	free(symbol);
#elif !defined(__CYGWIN__) && !defined(EMSCRIPTEN)
	size_t frames = backtrace(stack, max_stack_size);
	char** messages = backtrace_symbols(stack, frames);
	for (unsigned int i = 0; i < frames; i++) {
		std::string symbol_name = messages[i];
		if (symbol_name.find("get_stack_trace") != symbol_name.npos) {
			continue;
		}
		if (symbol_name.find("verbose_error") != symbol_name.npos) {
			continue;
		}

		auto begin_paren = symbol_name.find('(');
		auto plus = symbol_name.find('+');
		if (begin_paren != symbol_name.npos && plus != symbol_name.npos) {
			auto binary = symbol_name.substr(0, begin_paren);
			auto function = symbol_name.substr(begin_paren + 1, plus - begin_paren - 1);
			auto offset = symbol_name.substr(plus);

			char funcname[1024];
			std::size_t funcname_size = 1024;
			int status;
			char* demangled_name =
			    abi::__cxa_demangle(function.c_str(), funcname, &funcname_size, &status);

			sout += to_string(frames - i - 1, ": ", binary, " ");
			if (status == 0) {
				sout += demangled_name;
			} else {
				sout += function;
			}
			sout += to_string(" (", offset, "\n");
		} else {
			sout += to_string(frames - i - 1, ": ", symbol_name, "\n");
		}
	}
	std::free(messages);
#else
	sout += "<not available>\n";
#endif
	return sout;
}

// Removes the path from __FILE__ constants and keeps the name only.
std::string extract_file_name(const char* full_file_cstr) {
	using namespace std;

	// Extract the file name only.
	string file(full_file_cstr);
	auto pos = file.find_last_of("/\\");
	if (pos == string::npos) {
		pos = 0;
	}
	file = file.substr(pos + 1);  // Returns empty string if pos + 1 == length.

	return file;
}

void verbose_error_internal(const char* expr,
                            const char* function_name,
                            const char* full_file_cstr,
                            int line,
                            const char* args) {
	std::stringstream stream;
	stream << "Assumption failed: " << expr << " in " << function_name << " ("
	       << extract_file_name(full_file_cstr) << ":" << line << "). " << args;
	stream << get_stack_trace();
	throw std::logic_error(stream.str());
}
}  // namespace core
}  // namespace minimum
