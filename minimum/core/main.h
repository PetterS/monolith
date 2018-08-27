// Defines a main runner that captures exceptions, parses flags etc.
// The program should instead define this function:
//   int main_program(int num_args, char* args[]);
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

or simply

int main(int num_args, char* args[]) {
  return main_runner(main_program, num_args, args);
}
*/

#include <minimum/core/export.h>

MINIMUM_CORE_API int main_runner(int (*main_program)(int, char* []), int num_args, char* args[]);
