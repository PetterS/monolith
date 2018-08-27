#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <minimum/core/read_matrix.h>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *raw_data, size_t size) {
	std::string data((char *)raw_data, size);
	std::stringstream sin(data);
	try {
		auto mati = minimum::core::read_matrix<int>(sin, false);
		auto matd = minimum::core::read_matrix<double>(sin, false);
	} catch (...) {
		// Exceptions are OK.
	}
	return 0;
}
