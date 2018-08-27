#include <iomanip>
#include <iostream>

#include <minimum/ai/compute_equilibrium_sequential.h>
#include <minimum/ai/imperfect_games/goofspiel.h>
using namespace minimum::ai;
using namespace std;

void main_program() {
	typedef goofspiel::State State;

	auto result = equilibrium_sequential::compute<State>();
	equilibrium_sequential::print_result(result);
}

int main() {
	try {
		main_program();
		return 0;
	} catch (std::exception& error) {
		std::cerr << "ERROR: " << error.what() << std::endl;
		return 1;
	}
}
