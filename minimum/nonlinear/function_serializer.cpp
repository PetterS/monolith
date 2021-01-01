// Petter Strandmark.
#include <map>
#include <stdexcept>
#include <vector>

#include <minimum/nonlinear/function_serializer.h>

using namespace std;

namespace minimum {
namespace nonlinear {

std::ostream& operator<<(std::ostream& out, const Serialize& serializer) {
	if (!serializer.readonly_function) {
		throw runtime_error("Serializer << : Invalid function pointer.");
	}
	const Function& f = *serializer.readonly_function;
	f.write_to_stream(out);
	return out;
}

std::istream& operator>>(std::istream& in, const Serialize& serializer) {
	minimum_core_assert(serializer.writable_function, "Serializer << : Invalid function pointer.");
	minimum_core_assert(serializer.user_space, "Serializer << : Invalid user space pointer.");
	minimum_core_assert(serializer.factory, "Serializer << : Invalid factory pointer.");
	Function& f = *serializer.writable_function;
	f.read_from_stream(in, serializer.user_space, *serializer.factory);
	return in;
}
}  // namespace nonlinear
}  // namespace minimum
