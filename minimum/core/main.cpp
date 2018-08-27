// Defines a main runner that captures exceptions, parses flags etc.
// The program should instead define this function:
int main_program(int num_args, char* args[]);
// And a main function somewhat like this:
/*
int main(int num_args, char* args[]) {
  if (num_args >= 2 && minimum::core::ends_with(to_string(args[1]), "-help")) {
    gflags::SetUsageMessage(
      std::string("Runs colgen method for shift scheduling.\n   ")
      + std::string(args[0]) + string(" <input_file>"));
    gflags::ShowUsageWithFlagsRestrict(args[0], "colgen");
    return 0;
  }
  return main_runner(main_program, num_args, args);
}
*/

#include <iostream>

#ifndef EMSCRIPTEN
#include <absl/debugging/failure_signal_handler.h>
#include <absl/debugging/symbolize.h>
#endif

#include <gflags/gflags.h>

#include <minimum/core/flamegraph.h>
#include <minimum/core/main.h>
#include <minimum/core/scope_guard.h>
#include <minimum/core/string.h>
#include <version/version.h>

using namespace std;

DEFINE_string(flamegraph,
              "",
              "If set, prints flamegraph data to this filename. If set to \"stdout\" or \"-\", "
              "prints to console.");

int main_runner(int (*main_program)(int, char* []), int num_args, char* args[]) {
#ifndef EMSCRIPTEN
	absl::InitializeSymbolizer(args[0]);
	absl::FailureSignalHandlerOptions options;
	absl::InstallFailureSignalHandler(options);
#endif

	gflags::FlagSaver saver;  // Because Emscripten calls main multiple times.
	gflags::SetVersionString(version::revision);
	gflags::ParseCommandLineFlags(&num_args, &args, true);
	try {
		if (!FLAGS_flamegraph.empty()) {
			flamegraph::enable();
		}
		int code = main_program(num_args, args);
		// Only print flamegraph data if there was no exception.
		if (FLAGS_flamegraph == "stdout" || FLAGS_flamegraph == "-") {
			flamegraph::write_pretty(cout);
		} else if (!FLAGS_flamegraph.empty()) {
			flamegraph::write_to_file(FLAGS_flamegraph);
		}
		return code;
	} catch (exception& err) {
		cerr << "Error: " << err.what() << endl;
		return 1;
	} catch (...) {
		cerr << "Unknown exception." << endl;
		return 2;
	}
}
