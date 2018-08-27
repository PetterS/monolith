// E.g.
//
// extract_program output.umz.tmp program.umz 195 -1

#include <fstream>
#include <iostream>
using namespace std;

#include <minimum/core/string.h>
using namespace minimum::core;

int main(int argc, char* argv[]) {
	if (argc <= 4) {
		cerr << "Usage: " << argv[0]
		     << " <input data> <output program> <start offset> <end_offset>\n";
		return 100;
	}

	ifstream input_program(argv[1], ios::binary);
	ofstream binary_output(argv[2], ios::binary);
	int start_offset = from_string<int>(argv[3]);
	int end_offset = from_string<int>(argv[4]);
	if (!input_program) {
		cerr << "Could not open file.\n";
		return 1;
	}

	input_program.seekg(0, ios::end);
	auto file_size = input_program.tellg();
	if (end_offset < 0) {
		end_offset = file_size;
	}
	auto to_read = end_offset - start_offset;
	vector<char> data(to_read);
	input_program.seekg(start_offset);
	input_program.read(data.data(), to_read);
	if (!input_program) {
		cerr << "Read failed.\n";
		return 1;
	}
	binary_output.write(data.data(), to_read);
	if (!binary_output) {
		cerr << "Write failed.\n";
		return 2;
	}
	return 0;
}
