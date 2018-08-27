
#include <cstdint>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/record_stream.h>

namespace minimum {
namespace core {
void write_record(ostream* out, const std::string& data) {
	check(bool(*out), "Invalid stream.");

	if (data.size() < 127) {
		int8_t small_size = int8_t(data.size());
		out->write(reinterpret_cast<const char*>(&small_size), 1);
	} else {
		int8_t small_size = -1;
		int64_t big_size = int64_t(data.size());
		out->write(reinterpret_cast<const char*>(&small_size), 1);
		out->write(reinterpret_cast<const char*>(&big_size), sizeof(int64_t));
	}
	check(bool(*out), "Could not write size.");
	out->write(data.data(), data.size());
	check(bool(*out), "Could not write data.");
}

void read_record(istream* in, std::string* data) {
	check(bool(*in), "Invalid stream.");

	int64_t data_size = 0;
	int8_t small_size = -1;
	in->read(reinterpret_cast<char*>(&small_size), 1);
	if (small_size >= 0) {
		data_size = small_size;
	} else {
		check(small_size == -1, "Invalid small_size in file.");
		in->read(reinterpret_cast<char*>(&data_size), sizeof(int64_t));
	}
	check(bool(*in), "Could not read size.");
	check(data_size >= 0, "Invalid size in file.");
	data->resize(data_size);
	in->read(&(*data)[0], data_size);
	check(bool(*in), "Could not read data.");
}
}  // namespace core
}  // namespace minimum
