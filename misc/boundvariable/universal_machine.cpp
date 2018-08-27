// Decryption key for the Codex is (\b.bb)(\v.vv)06FHPVboundvarHRAk
//
// Start codex:
//    universal_machine codex.umz codex.out

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/time.h>
using namespace minimum::core;

#include <misc/boundvariable/universal_machine.h>

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " <input program> <binary output file>\n";
		return 1;
	}

	ifstream input_program(argv[1], ios::binary);
	ofstream binary_output(argv[2], ios::binary);
	try {
		Machine machine(input_program);
		machine.set_stdin(&cin);
		machine.set_stdout(&cout);
		machine.set_binary_out(&binary_output);
		machine.execute();
	} catch (std::exception& e) {
		cerr << "==> Error: " << e.what() << "\n";
	}
}