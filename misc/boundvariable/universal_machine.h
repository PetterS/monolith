#pragma once

#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

#include <minimum/core/check.h>
#include <minimum/core/time.h>
using namespace minimum::core;

struct Platter {
	uint32_t number;

	Platter(uint8_t A, uint8_t B, uint8_t C, uint8_t D) {
		number = (A << 24) + (B << 16) + (C << 8) + D;
	}

	Platter(uint32_t num) { number = num; }

	int operator_number() const { return (0b11110000000000000000000000000000 & number) >> 28; }

	int A() const { return (0b111000000 & number) >> 6; }

	int B() const { return (0b111000 & number) >> 3; }

	int C() const { return 0b111 & number; }

	int special_A() const { return (0b00001110000000000000000000000000 & number) >> 25; }

	uint32_t special_value() const { return 0b00000001111111111111111111111111 & number; }
};
static_assert(sizeof(Platter) == 4, "Wrong platter size.");

constexpr uint64_t LOWER_32 = 0xFFFFFFFFull;

class Machine {
   public:
	long long instructions_executed = 0;

	Machine(istream& input) : registers(8, 0) {
		original_platters.emplace_back();

		while (true) {
			uint8_t A, B, C, D;
			input.read((char*)&A, 1);
			input.read((char*)&B, 1);
			input.read((char*)&C, 1);
			input.read((char*)&D, 1);
			if (!input) {
				break;
			}
			original_platters[0].emplace_back(A, B, C, D);
		}
		cerr << "==> Read program of size " << original_platters[0].size() << ".\n";

		opcode_to_instruction.resize(100, nullptr);
		opcode_to_instruction[0] = &Machine::conditional_move;
		opcode_to_instruction[1] = &Machine::array_index;
		opcode_to_instruction[2] = &Machine::array_amendment;
		opcode_to_instruction[3] = &Machine::add;
		opcode_to_instruction[4] = &Machine::multiply;
		opcode_to_instruction[5] = &Machine::divide;
		opcode_to_instruction[6] = &Machine::not_and;
		opcode_to_instruction[7] = &Machine::halt;
		opcode_to_instruction[8] = &Machine::allocate;
		opcode_to_instruction[9] = &Machine::abandon;
		opcode_to_instruction[10] = &Machine::print;
		opcode_to_instruction[11] = &Machine::read;
		opcode_to_instruction[12] = &Machine::load_program;
		opcode_to_instruction[13] = &Machine::orthography;
	}

	Machine(const Machine& other) { *this = other; }

	Machine& operator=(const Machine& other) {
		// public members.
		instructions_executed = other.instructions_executed;
		// protected members.
		binary_output = other.binary_output;
		std_out = other.std_out;
		std_in = other.std_in;
		opcode_to_instruction = other.opcode_to_instruction;

		platters = other.platters;
		registers = other.registers;
		free_platters = other.free_platters;
		finger = other.finger;

		// Compile to fill in the other state.
		compile();
		return *this;
	}

	virtual ~Machine() {}

	void set_stdin(istream* in) { std_in = in; }

	void set_stdout(ostream* out) { std_out = out; }

	void set_binary_out(ostream* in) { binary_output = in; }

	void clear() {
		platters = original_platters;
		compile();

		finger = 0;
		for (auto& reg : registers) {
			reg = 0;
		}
		use_binary_output = false;
		instructions_executed = 0;
	}

	class HaltException {};

	void execute() {
		clear();
		start_time = wall_time();
		resume();
		double elapsed_time = wall_time() - start_time;
		cerr << "\n==> Finished " << instructions_executed << " instructions in " << elapsed_time
		     << " seconds.\n";
	}

	void resume() {
		try {
			while (true) {
				// Since this is an extremely hot loop, the universal machine
				// can be sped up by removing the safety checks here.
				// instructions_executed++;
				// check(finger < instructions.size(), "Invalid finger position.");
				Instruction instruction = instructions[finger];
				// check(instruction != nullptr, "Invalid instruction.");

				// Run the current instruction.
				(this->*instruction)();
			}
		} catch (HaltException&) {
			// Done.
		}
	}

   protected:
	ostream* binary_output = nullptr;
	ostream* std_out = nullptr;
	istream* std_in = nullptr;

	using Instruction = void (Machine::*)();
	vector<Instruction> opcode_to_instruction;

	// All memory of the machine.
	vector<vector<Platter>> platters;
	vector<uint32_t> registers;
	// Indices of free arrays that can be reused.
	vector<uint32_t> free_platters;
	// Program position.
	uint32_t finger = 0;

	// Compiled array 0.
	vector<uint32_t*> registerA;
	vector<uint32_t*> registerB;
	vector<uint32_t*> registerC;
	vector<Instruction> instructions;
	vector<uint32_t> special_value;

	void conditional_move() {
		if (*registerC[finger] != 0) {
			*registerA[finger] = *registerB[finger];
		}
		finger++;
	}

	void array_index() {
		auto array_number = *registerB[finger];
		auto offset = *registerC[finger];
		minimum_core_assert(array_number < platters.size());
		minimum_core_assert(offset < platters[array_number].size());
		*registerA[finger] = platters[array_number][offset].number;
		finger++;
	}

	void array_amendment() {
		auto array_number = *registerA[finger];
		auto offset = *registerB[finger];
		minimum_core_assert(array_number < platters.size());
		minimum_core_assert(offset < platters[array_number].size());
		platters[array_number][offset].number = *registerC[finger];
		if (array_number == 0) {
			compile(offset);
		}
		finger++;
	}

	void add() {
		uint64_t result = uint64_t(*registerB[finger]) + uint64_t(*registerC[finger]);
		*registerA[finger] = LOWER_32 & result;
		finger++;
	}

	void multiply() {
		uint64_t result = uint64_t(*registerB[finger]) * uint64_t(*registerC[finger]);
		*registerA[finger] = LOWER_32 & result;
		finger++;
	}

	void divide() {
		minimum_core_assert(*registerC[finger] != 0);
		*registerA[finger] = *registerB[finger] / *registerC[finger];
		finger++;
	}

	void not_and() {
		*registerA[finger] = ~(*registerB[finger] & *registerC[finger]);
		finger++;
	}

	void halt() { throw HaltException(); }

	void allocate() {
		auto capacity = *registerC[finger];
		if (!free_platters.empty()) {
			*registerB[finger] = free_platters.back();
			free_platters.pop_back();
			platters[*registerB[finger]].resize(capacity, Platter(0));
		} else {
			platters.emplace_back(capacity, Platter(0));
			*registerB[finger] = platters.size() - 1;
		}
		finger++;
	}

	void abandon() {
		auto array_number = *registerC[finger];
		minimum_core_assert(array_number < platters.size());
		platters[array_number].clear();
		free_platters.push_back(array_number);
		finger++;
	}

	void print() {
		minimum_core_assert(0 <= *registerC[finger] && *registerC[finger] < 256);
		uint8_t byte = *registerC[finger];
		if (binary_output != nullptr) {
			binary_output->write((char*)&byte, 1);
		}
		if (!use_binary_output && ((byte < 0x20 && byte != '\n' && byte != '\t') || byte >= 128)) {
			cerr << "\n==> Non-printable ASCII found. Supressing output to stdout until program "
			        "halts.\n";
			use_binary_output = true;
		}
		if (!use_binary_output && std_out != nullptr) {
			*std_out << byte;
		}
		finger++;
	}

	void read() {
		if (std_in == nullptr || !*std_in) {
			*registerC[finger] = ~(uint32_t)0;
		} else {
			int ch = std_in->get();
			minimum_core_assert(0 <= ch && ch < 256);
			if (ch == '`') {
				dump();
			}
			*registerC[finger] = ch;
		}
		finger++;
	}

	void load_program() {
		auto array_number = *registerB[finger];
		finger = *registerC[finger];
		if (array_number != 0) {
			minimum_core_assert(array_number < platters.size());
			platters[0] = platters[array_number];
			compile();
		}
	}

	void orthography() {
		*registerA[finger] = special_value[finger];
		finger++;
	}

	void dump() {
		int num_platters = 0;
		int platters_allocated = 0;
		for (auto& array : platters) {
			platters_allocated += array.size();
			if (array.size() > 0) {
				num_platters++;
			}
		}
		cerr << "\n";
		cerr << "==> DEBUG\n";
		cerr << "==> " << num_platters << " arrays exist with " << platters_allocated
		     << " platters.\n";
		cerr << "==> " << instructions_executed << " instructions executed.\n";
	}

   private:
	void compile() {
		auto size = platters[0].size();
		registerA.resize(size);
		registerB.resize(size);
		registerC.resize(size);
		instructions.resize(size);
		special_value.resize(size, 0);
		for (size_t p = 0; p < size; ++p) {
			compile(p);
		}
	}

	void compile(size_t p) {
		auto op = platters[0][p].operator_number();
		registerA[p] = &registers[platters[0][p].A()];
		registerB[p] = &registers[platters[0][p].B()];
		registerC[p] = &registers[platters[0][p].C()];
		instructions[p] = opcode_to_instruction[op];

		if (op == 13) {
			registerA[p] = &registers[platters[0][p].special_A()];
			special_value[p] = platters[0][p].special_value();
		}
	}

	bool use_binary_output = false;

	// Statistics.
	double start_time = 0;

	vector<vector<Platter>> original_platters;
};
