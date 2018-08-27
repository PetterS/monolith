#include <iomanip>
#include <iostream>

#include <minimum/ai/compute_equilibrium_sequential.h>
#include <minimum/ai/imperfect_games/leduc_holdem.h>
using namespace minimum::ai;
using namespace std;

void main_program() {
	typedef leduc_holdem::State State;
	typedef leduc_holdem::InformationSet InformationSet;

	auto result = equilibrium_sequential::compute<State>();
	if (false) {
		// Print all strategies.
		equilibrium_sequential::print_result(result);
	}

	cout << "Game value is " << setprecision(14) << result.value << " for player 0.\n\n";

	cout << "FIRST PLAYER: First move.\n";
	for (int card = 1; card <= 3; ++card) {
		cout << "When holding card " << card << "\n";

		State state(card, (card % 3) + 1, card);
		auto policy = result.player_strategies[0].at(InformationSet(state));
		for (auto move : policy) {
			cout << "-- " << State::move_name(move.first) << " with prob. " << move.second << "\n";
		}
	}
	cout << "\n";

	cout << "SECOND PLAYER: Response after first player raises.\n";
	for (int card = 1; card <= 3; ++card) {
		cout << "When holding card " << card << "\n";

		State state((card % 3) + 1, card, card);
		state.move(State::Move::raise);

		auto policy = result.player_strategies[1].at(InformationSet(state));
		for (auto move : policy) {
			cout << "-- " << State::move_name(move.first) << " with prob. " << move.second << "\n";
		}
	}
	cout << "\n";

	cout << "FIRST PLAYER: Response after second player re-raises.\n";
	for (int card = 1; card <= 3; ++card) {
		cout << "When holding card " << card << "\n";

		State state(card, (card % 3) + 1, card);
		state.move(State::Move::raise);
		state.move(State::Move::raise);

		auto policy = result.player_strategies[0].at(InformationSet(state));
		for (auto move : policy) {
			cout << "-- " << State::move_name(move.first) << " with prob. " << move.second << "\n";
		}
	}
	cout << "\n";

	cout << "FIRST PLAYER: Response after a first round of (raise, raise, call).\n";
	for (int card = 1; card <= 3; ++card) {
		for (int table_card = 1; table_card <= 3; ++table_card) {
			cout << "When holding card " << card << " and the community card is " << table_card
			     << "\n";

			State state(card, (card % 3) + 1, table_card);
			// First round.
			state.move(State::Move::raise);
			state.move(State::Move::raise);
			state.move(State::Move::call);

			auto policy = result.player_strategies[0].at(InformationSet(state));
			for (auto move : policy) {
				cout << "-- " << State::move_name(move.first) << " with prob. " << move.second
				     << "\n";
			}
		}
	}
	cout << "\n";

	cout
	    << "SECOND PLAYER: Response after a first round of (raise, raise, call), then first player "
	       "raises.\n";
	for (int card = 1; card <= 3; ++card) {
		for (int table_card = 1; table_card <= 3; ++table_card) {
			cout << "When holding card " << card << " and the community card is " << table_card
			     << "\n";

			State state((card % 3) + 1, card, table_card);
			// First round.
			state.move(State::Move::raise);
			state.move(State::Move::raise);
			state.move(State::Move::call);
			// Second round.
			state.move(State::Move::raise);

			auto policy = result.player_strategies[1].at(InformationSet(state));
			for (auto move : policy) {
				cout << "-- " << State::move_name(move.first) << " with prob. " << move.second
				     << "\n";
			}
		}
	}
	cout << "\n";

	// The following policy for player 0 looks incorrect:
	//
	// -- At information set Round 1, 2 raises, P0: $8, P1 : $12.Card : 3, flop : 3. History : call
	// raise raise call raise raise :
	// ---- Do fold with probability 0.113846
	// ---- Do call with probability 0.886154
	//
	// Because (3, 3) is a guaranteed win. But the equilibrium contains the following entries for
	// player 1:
	//
	// -- At information set Round 1, 1 raises, P0: $8, P1 : $4.Card : 1, flop : 3. History : call
	// raise raise call raise :
	// ---- Do fold with probability 1
	// ---- Do raise with probability 0
	// ---- Do call with probability 0
	// -- At information set Round 1, 1 raises, P0: $8, P1 : $4.Card : 2, flop : 3. History : call
	// raise raise call raise :
	// ---- Do call with probability 0.227203
	// ---- Do raise with probability 0
	// ---- Do fold with probability 0.772797
	//
	// So (call raise raise call raise raise) will in fact never be played. And player 1 will be
	// worse off if she starts raising on (call raise raise call raise), since player 0 will still
	// call the raise with probability ≈0.9.
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
