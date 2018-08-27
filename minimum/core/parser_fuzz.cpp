#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <string_view>
#include <vector>

#include <minimum/core/parser.h>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *raw_data, size_t size) {
	std::string_view expr((const char *)raw_data, size);
	try {
		minimum::core::Expression expression = minimum::core::Parser(expr).parse();
		expression.evaluate<double>({});
	} catch (std::runtime_error &) {
		// Exceptions are OK.
	}
	return 0;
}
