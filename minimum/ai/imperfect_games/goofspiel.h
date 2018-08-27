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

namespace goofspiel {

constexpr int num_cards = 4;

class InformationSet;

class State {
   public:
	friend class InformationSet;
	typedef goofspiel::InformationSet InformationSet;

	typedef int Move;
	typedef std::array<double, 2> Reward;

	static int num_players() { return 2; }

	static std::string move_name(Move move) { return minimum::core::to_string(move); }

	// Returns a vector of all possible initial states (with equal probability).
	// This is needed for solving the game exactly.
	static std::vector<State> all_initial_states() {
		std::vector<State> states(1);
		return states;
	}

	State() {
		for (auto card : range(1, num_cards + 1)) {
			player_cards[0].push_back(card);
			player_cards[1].push_back(card);
		}
	}

	// Whether this is a terminal state where the game has ended.
	bool terminal() const { return player_cards[1].empty(); }

	// The player that should perform a move.
	int player() const { return player_to_move; }

	// Computes which player is the winner and returns the signed reward
	// for both players.
	Reward reward() const {
		minimum_core_assert(terminal());
		return {game_result, -game_result};
	}

	// All moves that are possible in the state.
	std::vector<Move> possible_moves() const {
		minimum_core_assert(!terminal());
		return player_cards[player_to_move];
	}

	// Performs a move for the current player and updates the state in-place
	// to the resulting new state.
	State& move(Move move) {
		minimum_core_assert(!terminal());
		played_card[player_to_move] = move;

		if (player_to_move == 1) {
			double card_value = num_cards - player_cards[0].size() + 1;

			for (int i : range(2)) {
				auto itr =
				    std::find(player_cards[i].begin(), player_cards[i].end(), played_card[i]);
				minimum_core_assert(itr != player_cards[i].end());
				player_cards[i].erase(itr);
				minimum_core_assert(played_card[i] > 0);
				history[i].push_back(played_card[i]);
			}

			if (played_card[0] > played_card[1]) {
				game_result += card_value;
			} else if (played_card[0] < played_card[1]) {
				game_result -= card_value;
			} else {
				// Players split the score.
			}

			minimum_core_assert(player_cards[0].size() == player_cards[1].size());

			// Play the last card.
			if (player_cards[0].size() == 1) {
				if (player_cards[0].back() > player_cards[1].back()) {
					game_result += num_cards;
				} else if (player_cards[0].back() < player_cards[1].back()) {
					game_result -= num_cards;
				} else {
					// Players split the score.
				}
				player_cards[0].clear();
				player_cards[1].clear();
			}
		}
		player_to_move = 1 - player_to_move;
		return *this;
	}

   private:
	std::vector<int> player_cards[2];
	int played_card[2] = {-1, -1};
	std::vector<int> history[2];
	int player_to_move = 0;
	double game_result = 0;
};

class InformationSet {
   public:
	// Creates the "null" information sets. Can not be equal to any other.
	InformationSet() {
		history[0].push_back(-1);
		player_cards[0].push_back(-1);
	}

	// Creates an InformationSet from a State, using all the information the player
	// has access to.
	InformationSet(const State& state) {
		player_cards[0] = state.player_cards[0];
		player_cards[1] = state.player_cards[1];
		history[0] = state.history[0];
		history[1] = state.history[1];
	}

	// Used when storing information sets as keys in hash maps.
	bool operator==(const InformationSet& other) const {
		return history[0] == other.history[0] && history[1] == other.history[1];

		// The error below makes it a game of imperfect recall.
		// return player_cards[0] == other.player_cards[0] && player_cards[1] ==
		// other.player_cards[1];
	}

	// String respresentation for writing down strategies.
	std::string str() const {
		std::stringstream sout;
		for (int i : range(2)) {
			sout << "Player " << i << " has played " << to_string(history[i]) << ".";
		}
		return sout.str();
	}

	std::vector<int> player_cards[2];
	std::vector<int> history[2];
};
}  // namespace goofspiel

namespace std {
template <>
struct hash<goofspiel::InformationSet> {
	std::hash<int> hasher;
	size_t operator()(goofspiel::InformationSet const& is) const {
		size_t h = 0;
		for (int i : range(2)) {
			for (auto card : is.history[i]) {
				h = hash_combine(h, hasher(card));
			}
		}
		return h;
	}
};
}  // namespace std
