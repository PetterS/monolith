
#include <iostream>
#include <unordered_map>

#include <minimum/ai/compute_equilibrium.h>
#include <minimum/ai/compute_equilibrium_sequential.h>
#include <minimum/ai/imperfect_games/kuhn_poker.h>
#include <minimum/ai/imperfect_games/rock_paper_scissors.h>
#include <minimum/ai/imperfect_games/stengel.h>
using namespace minimum::ai;
using namespace std;

template <typename State>
void print_strategies() {
	auto result = equilibrium_sequential::compute<State>();
	equilibrium_sequential::print_result(result);
}

int main() {
	try {
		cout << "ROCK, PAPER, SCISSORS\n";
		print_strategies<rock_paper_scissors::State>();
		cout << "ROCK, PAPER, SCISSORS, SPOCK, LIZARD\n";
		print_strategies<rock_paper_scissors::SpockState>();
		cout << "\nKUHN POKER\n";
		print_strategies<kuhn_poker::State<>>();
		cout << "\nSTENGEL\n";
		print_strategies<stengel::State>();
	} catch (std::exception& error) {
		std::cerr << "ERROR: " << error.what() << std::endl;
		return 1;
	}
}
