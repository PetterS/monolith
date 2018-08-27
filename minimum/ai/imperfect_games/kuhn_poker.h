#pragma once

#include <array>
#include <functional>
#include <random>
#include <sstream>

#include <minimum/core/check.h>
#include <minimum/core/hash.h>

namespace kuhn_poker {

template <typename State>
class InformationSet;

enum class Move : signed char { null = 0, bet, check, call, fold };

template <int number_of_cards = 3>
class State {
   public:
	friend class InformationSet<State>;
	typedef kuhn_poker::InformationSet<State> InformationSet;

	typedef kuhn_poker::Move Move;
	typedef std::array<double, 2> Reward;

	static const int ante_size = 1;
	static const int bet_size = 1;

	static int num_players() { return 2; }

	static std::string move_name(Move move) {
		if (move == Move::null) {
			return "null";
		} else if (move == Move::bet) {
			return "bet";
		} else if (move == Move::check) {
			return "check";
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
		for (int c1 = 1; c1 <= number_of_cards; ++c1) {
			for (int c2 = 1; c2 <= number_of_cards; ++c2) {
				if (c1 != c2) {
					states.emplace_back(c1, c2);
				}
			}
		}
		return states;
	}

	template <typename Engine>
	State(Engine* random_engine) {
		std::uniform_int_distribution<int> cards(1, number_of_cards);
		player0_card = cards(*random_engine);
		do {
			player1_card = cards(*random_engine);
		} while (player0_card == player1_card);
	}

	State(int player0_card, int player1_card) {
		minimum_core_assert(player0_card != player1_card);
		minimum_core_assert(1 <= player0_card && player0_card <= number_of_cards);
		minimum_core_assert(1 <= player1_card && player1_card <= number_of_cards);
		this->player0_card = player0_card;
		this->player1_card = player1_card;
	}

	bool terminal() const { return game_result != 0; }

	int player() const { return game_round % 2; }

	Reward reward() const {
		minimum_core_assert(terminal());
		return {double(game_result), double(-game_result)};
	}

	std::vector<Move> possible_moves() const {
		minimum_core_assert(!terminal());

		if (game_round == 0) {
			return {Move::check, Move::bet};
		} else if (game_round == 1) {
			if (first_move == Move::check) {
				return {Move::check, Move::bet};
			} else {
				return {Move::fold, Move::call};
			}
		} else if (game_round == 2) {
			return {Move::fold, Move::call};
		} else {
			minimum_core_assert(false && "Incorrect game round.");
			return {};
		}
	}

	State& move(Move move) {
		minimum_core_assert(!terminal());

		// Too much overhead outside testing.
		// auto moves = possible_moves();
		// minimum_core_assert(find(begin(moves), end(moves), move) != end(moves));

		if (game_round == 0) {
			// Player 0 can check or bet 1.
			minimum_core_assert(move == Move::check || move == Move::bet);
			first_move = move;
		} else if (game_round == 1) {
			if (first_move == Move::check) {
				// If player 0 checks then player 2 can check or bet 1.
				minimum_core_assert(move == Move::check || move == Move::bet);
				if (move == Move::check) {
					showdown(ante_size);
				}
			} else {
				// If player 0 bets then player 2 can fold or call.
				minimum_core_assert(move == Move::fold || move == Move::call);
				if (move == Move::fold) {
					game_result = ante_size;
				} else {
					showdown(ante_size + bet_size);
				}
			}
		} else if (game_round == 2) {
			// If player 0 checks and player 2 bets, player 0 can fold or call.
			minimum_core_assert(first_move == Move::check);
			minimum_core_assert(move == Move::fold || move == Move::call);
			if (move == Move::fold) {
				game_result = -ante_size;
			} else {
				showdown(ante_size + bet_size);
			}
		}

		game_round++;
		return *this;
	}

   private:
	void showdown(int player_bet) {
		if (player0_card > player1_card) {
			game_result = player_bet;
		} else {
			game_result = -player_bet;
		}
	}

	signed char player0_card = -1;
	signed char player1_card = -1;

	signed char game_round = 0;
	signed char game_result = 0;
	Move first_move = Move::null;
};

template <typename State>
class InformationSet {
   public:
	typedef typename State::Move Move;

	InformationSet() : game_round(-1), my_card(-1), this_player(-1) {}

	InformationSet(const State& state)
	    : game_round(state.game_round), this_player(state.player()), first_move(state.first_move) {
		if (player() == 0) {
			my_card = state.player0_card;
		} else {
			my_card = state.player1_card;
		}
	}

	int player() const { return this_player; }
	int card() const { return my_card; }
	int round() const { return game_round; }

	bool operator==(const InformationSet& other) const {
		return this_player == other.this_player && my_card == other.my_card
		       && game_round == other.game_round && first_move == other.first_move;
	}

	std::string str() const {
		std::stringstream sout;
		sout << "Player:" << player() << ", card:" << card() << " round:" << round();
		if (first_move != Move::null) {
			sout << " 1stmove:" << State::move_name(first_move);
		}
		return sout.str();
	}

	signed char game_round;
	signed char my_card;
	signed char this_player;
	Move first_move = Move::null;
};

}  // namespace kuhn_poker

namespace std {
template <int number_of_cards>
struct hash<kuhn_poker::InformationSet<kuhn_poker::State<number_of_cards>>> {
	size_t operator()(
	    const kuhn_poker::InformationSet<kuhn_poker::State<number_of_cards>>& is) const {
		return minimum::core::hasher(
		    is.this_player, is.my_card, is.game_round, (signed char)is.first_move);
	}
};

template <>
struct hash<typename kuhn_poker::Move> {
	std::hash<signed char> hasher;
	size_t operator()(const kuhn_poker::Move& move) const { return hasher((unsigned char)move); }
};
}  // namespace std
