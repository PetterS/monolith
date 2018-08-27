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

class BreakpointException {};

class BreakerMachine : public Machine {
   public:
	BreakerMachine(istream& input_program) : Machine(input_program) {
		opcode_to_instruction[10] = static_cast<Instruction>(&BreakerMachine::my_print);
		opcode_to_instruction[11] = static_cast<Instruction>(&BreakerMachine::my_read);
	}

	bool trap = false;
	bool login_ok = false;
	long instructions_until_trap = 0;

	void set_trap() {
		trap = true;
		login_ok = false;
		instructions_until_trap = 0;
	}

	void my_print() {
		if (trap) {
			if (*registerC[finger] == 'y') {
				login_ok = true;
				// throw BreakpointException();
			} else if (*registerC[finger] == 'n') {
				login_ok = false;
				// throw BreakpointException();
			}
		}
		Machine::print();
	}

	void my_read() {
		if (std_in != nullptr && std_in->peek() == '`') {
			throw BreakpointException();
		}
		Machine::read();
	}
};

bool try_pass(BreakerMachine& breaker, const string& user, const string& password) {
	istringstream cmd("./trypass.exe " + user + " " + password + "\n`", ios::binary);
	breaker.set_trap();
	breaker.login_ok = false;
	breaker.set_stdin(&cmd);
	try {
		breaker.resume();
	} catch (BreakpointException& e) {
	}
	return breaker.login_ok;
}

auto words = {"airplane",
              "alphabet",
              "aviator",
              "bidirectional",
              "changeme",
              "creosote",
              "cyclone",
              "december",
              "dolphin",
              "elephant",
              "ersatz",
              "falderal",
              "functional",
              "future",
              "guitar",
              "gymnast",
              "hello",
              "imbroglio",
              "january",
              "joshua",
              "kernel",
              "kingfish",
              "(\\b.bb)(\\v.vv)",
              "millennium",
              "monday",
              "nemesis",
              "oatmeal",
              "october",
              "paladin",
              "pass",
              "password",
              "penguin",
              "polynomial",
              "popcorn",
              "qwerty",
              "sailor",
              "swordfish",
              "symmetry",
              "system",
              "tattoo",
              "thursday",
              "tinman",
              "topography",
              "unicorn",
              "vader",
              "vampire",
              "viper",
              "warez",
              "xanadu",
              "xyzzy",
              "zephyr",
              "zeppelin",
              "zxcvbnm"};

int main_program(int argc, char* argv[]) {
	if (argc < 3) {
		cerr << "Usage: " << argv[0] << " <operating system> <user>\n";
		return 1;
	}
	ifstream input_program(argv[1], ios::binary);
	BreakerMachine machine(input_program);
	machine.set_stdout(&cout);

	string startup = "guest\n";
	startup += "mail\n";
	startup += "/bin/umodem trypass.bas STOP\n";
	startup += "X IF CHECKPASS(ARG(II), ARG(III)) THEN GOTO XXXX\n";
	startup += "XX PRINT \"no\" + CHR(X)\n";
	startup += "XXX END\n";
	startup += "XXXX PRINT \"yes\" + CHR(X)\n";
	startup += "STOP\n";
	startup += "/bin/qbasic trypass.bas\n";
	startup += "`";

	try {
		istringstream startup_input(startup, ios::binary);
		machine.set_stdin(&startup_input);
		machine.execute();
	} catch (BreakpointException& e) {
	}
	machine.instructions_executed = 0;
	cerr << "\n==> Machine is in read state.\n";
	machine.set_stdout(nullptr);

	for (string user : {"ftd", "knr", "gardener", "ohmega", "yang", "howie", "hmonk", "bbarker"}) {
		for (auto word : words) {
			if (try_pass(machine, user, word)) {
				cerr << "==> Found " << user << ":" << word << "\n";
			}

			for (char ch1 = '0'; ch1 <= '9'; ++ch1) {
				for (char ch2 = '0'; ch2 <= '9'; ++ch2) {
					string password = word;
					password += ch1;
					password += ch2;
					if (try_pass(machine, user, password)) {
						cerr << "==> Found " << user << ":" << password << "\n";
					}
				}
			}
		}
	}

	return 1;
}

int main(int argc, char* argv[]) {
	try {
		return main_program(argc, argv);
	} catch (std::exception& e) {
		cerr << "==> Error: " << e.what() << "\n";
	}
}