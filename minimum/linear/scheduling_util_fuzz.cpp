#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>

#include <minimum/linear/scheduling_util.h>

using namespace std;

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* raw_data, size_t size) {
	string data((char*)raw_data, size);
	istringstream sin(data);
	try {
		auto problem = read_nottingham_instance(sin, false);
		// cout << "Number of vars: " << ip->get_number_of_variables() << "\n";
	} catch (std::logic_error&) {
		// From minimum_core_assert
	} catch (std::runtime_error&) {
		// From check
	}
	return 0;
}
