#include <iomanip>
#include <iostream>

#include <minimum/ai/compute_equilibrium_sequential.h>
#include <minimum/ai/imperfect_games/kuhn_poker.h>
using namespace minimum::ai;
using namespace std;

void main_program() {
	// 3 cards is the canonical Kuhn Poker.
	constexpr int num_cards = 52;

	typedef kuhn_poker::State<num_cards> State;
	typedef kuhn_poker::InformationSet<State> InformationSet;

	auto result = equilibrium_sequential::compute<State>();
	cout << "Game value is " << result.value << " for player 0.\n\n";

	for (int player0_card = 1; player0_card <= num_cards; ++player0_card) {
		State state(player0_card, player0_card % num_cards + 1);
		InformationSet is(state);

		auto& strategy = result.player_strategies[0].at(is);
		for (auto& move : strategy) {
			if (move.first == State::Move::bet) {
				std::cout << "With card " << std::setw(2) << player0_card
				          << ", bet with probability " << move.second << "\n";
			}
		}
	}
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
