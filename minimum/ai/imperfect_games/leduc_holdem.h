#pragma once

#include <array>
#include <functional>
#include <random>
#include <sstream>

#include <minimum/core/check.h>
#include <minimum/core/hash.h>
#include <minimum/core/range.h>
using minimum::core::hash_combine;
using minimum::core::range;

namespace leduc_holdem {

enum class Move : signed char { null = 0, raise, call, fold };
const std::vector<int> all_cards = {1, 1, 2, 2, 3, 3};
const int ante_size = 1;
const int bet_size[2] = {2, 4};

class InformationSet;

class State {
   public:
	friend class InformationSet;
	typedef leduc_holdem::InformationSet InformationSet;

	typedef leduc_holdem::Move Move;
	typedef std::array<double, 2> Reward;

	static int num_players() { return 2; }

	static std::string move_name(Move move) {
		if (move == Move::null) {
			return "null";
		} else if (move == Move::raise) {
			return "raise";
		} else if (move == Move::call) {
			return "call";
		} else if (move == Move::fold) {
			return "fold";
		} else {
			throw std::invalid_argument("Invalid move.");
		}
	}

	// Returns a vector of all possible initial states (with equal probability).
	// This is needed for solving the game exactly.
	static std::vector<State> all_initial_states() {
		std::vector<State> states;
		for (auto i : range(all_cards.size())) {
			for (auto j : range(all_cards.size())) {
				if (i == j) {
					continue;
				}
				for (auto k : range(all_cards.size())) {
					if (i == k || j == k) {
						continue;
					}
					states.emplace_back(all_cards[i], all_cards[j], all_cards[k]);
				}
			}
		}
		return states;
	}

	State(int player0_card_, int player1_card_, int flop_)
	    : player0_card(player0_card_), player1_card(player1_card_), flop(flop_) {}

	// Whether this is a terminal state where the game has ended.
	bool terminal() const { return game_round >= 2; }

	// The player that should perform a move.
	int player() const { return player_to_move; }

	// Computes which player is the winner and returns the signed reward
	// for both players.
	Reward reward() const {
		minimum_core_assert(terminal());
		return {double(game_result), double(-game_result)};
	}

	// All moves that are possible in the state.
	std::vector<Move> possible_moves() const {
		minimum_core_assert(!terminal());

		if (player_bets[0] == player_bets[1]) {
			return {Move::raise, Move::call};
		} else if (raises_in_round < 2) {
			return {Move::raise, Move::call, Move::fold};
		} else {
			return {Move::call, Move::fold};
		}
	}

	// Performs a move for the current player and updates the state in-place
	// to the resulting new state.
	State& move(Move move) {
		minimum_core_assert(!terminal());
		all_moves.emplace_back(move);
		// std::cout << "P" << int(player_to_move) << ", " << move_name(move) << "\n";

		if (move == Move::call) {
			minimum_core_assert(player_bets[player_to_move] <= player_bets[1 - player_to_move]);

			if (player_to_move == 0 && player_bets[0] == player_bets[1]) {
				// Nothing happens. It’s the next player’s turn.
			} else {
				// Put in the amount of money owed the pot.
				auto needed = std::max(player_bets[0], player_bets[1]);
				player_bets[0] = needed;
				player_bets[1] = needed;

				// Go to next round.
				game_round++;
				raises_in_round = 0;
				player_to_move = 1;  // Will then be set to 0 just before returning.

				if (terminal()) {
					showdown();
				}
			}
		} else if (move == Move::raise) {
			minimum_core_assert(raises_in_round < 2);
			raises_in_round++;

			// First, put in the amount of money owed the pot.
			auto needed = std::max(player_bets[0], player_bets[1]);
			if (needed == 1) {
				// If there is only the antes on the table, they should
				// be included in the first bet.
				needed = 0;
			}
			player_bets[player_to_move] = needed;

			// Then, put in an extra amount.
			player_bets[player_to_move] += bet_size[game_round];
		} else if (move == Move::fold) {
			minimum_core_assert(player_bets[player_to_move] < player_bets[1 - player_to_move]);
			if (player_to_move == 0) {
				game_result = -player_bets[0];
			} else {
				game_result = player_bets[1];
			}
			game_round = 2;
		} else {
			minimum_core_assert(false, "Invalid move.");
		}

		player_to_move = 1 - player_to_move;
		return *this;
	}

   private:
	// When both players call, the cards are shown and the
	// winner takes the pot.
	void showdown() {
		// At showdown, both players have put an equal amount.
		minimum_core_assert(player_bets[0] == player_bets[1]);

		if (player0_card == flop) {
			game_result = player_bets[1];
		} else if (player1_card == flop) {
			game_result = -player_bets[0];
		} else if (player0_card > player1_card) {
			game_result = player_bets[1];
		} else if (player1_card > player0_card) {
			game_result = -player_bets[0];
		} else {
			// Split pot.
			game_result = 0;
		}
	}

	const signed char player0_card;
	const signed char player1_card;
	const signed char flop;  // Card shown on the table in the second round.

	signed char game_round = 0;
	signed char raises_in_round = 0;  // Number of raises in the currect betting round.
	signed char player_to_move = 0;
	signed char player_bets[2] = {ante_size,
	                              ante_size};  // Total amount on the table for each player.

	std::vector<Move>
	    all_moves;  // Not used for game logic, but used for defining the information set.
	signed char game_result = 0;  // Updated during showdown or when a player folds.
};

class InformationSet {
   public:
	// Creates the "null" information sets. Can not be equal to any other.
	InformationSet() {}

	// Creates an InformationSet from a State, using all the information the player
	// has access to.
	InformationSet(const State& state)
	    : all_moves(state.all_moves),
	      game_round(state.game_round),
	      raises_in_round(state.raises_in_round),
	      player_to_move(state.player_to_move) {
		player_bets[0] = state.player_bets[0];
		player_bets[1] = state.player_bets[1];

		if (state.player_to_move == 0) {
			my_card = state.player0_card;
		} else {
			my_card = state.player1_card;
		}

		// In the first round, we do not know the flop, but in
		// the second, we do.
		if (state.game_round >= 1) {
			flop = state.flop;
		} else {
			flop = -1;
		}
	}

	// Used when storing information sets as keys in hash maps.
	bool operator==(const InformationSet& other) const {
		return all_moves == other.all_moves && my_card == other.my_card && flop == other.flop;
	}

	// String respresentation for writing down strategies.
	std::string str() const {
		std::stringstream sout;
		sout << "Round " << int(game_round) << ", ";
		sout << int(raises_in_round) << " raises, ";
		sout << "P0: $" << int(player_bets[0]) << ", ";
		sout << "P1: $" << int(player_bets[1]) << ". ";
		sout << "Card: " << int(my_card) << ", ";
		sout << "flop: " << int(flop) << ". ";
		sout << "History: ";
		for (auto move : all_moves) {
			sout << State::move_name(move) << " ";
		}
		return sout.str();
	}

	// These are used to define the information state for
	// the purposes of the game.
	const std::vector<Move> all_moves;
	signed char my_card = -1;
	signed char flop = -1;

	// These are only used to provide a better, human-readable
	// description string.
	signed char game_round = -1;
	signed char raises_in_round = -1;
	signed char player_to_move = -1;
	signed char player_bets[2] = {-1, -1};
};
}  // namespace leduc_holdem

namespace std {
template <>
struct hash<leduc_holdem::InformationSet> {
	std::hash<signed char> hasher;
	size_t operator()(leduc_holdem::InformationSet const& is) const {
		size_t h = hasher(is.my_card);
		h = hash_combine(h, hasher(is.flop));
		for (auto move : is.all_moves) {
			h = hash_combine(h, hasher((signed char)move));
		}
		return h;
	}
};

template <>
struct hash<leduc_holdem::Move> {
	std::hash<signed char> hasher;
	size_t operator()(leduc_holdem::Move const& move) const { return hasher((unsigned char)move); }
};
}  // namespace std
