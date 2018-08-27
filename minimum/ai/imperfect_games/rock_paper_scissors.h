#ifndef ROCK_PAPER_SCISSORS_HEADER
#define ROCK_PAPER_SCISSORS_HEADER

#include <array>
#include <functional>
#include <random>
#include <sstream>

#include <minimum/core/check.h>

namespace rock_paper_scissors {

template <typename State>
class InformationSet;

template <bool include_spock>
class TemplatedState {
   public:
	typedef rock_paper_scissors::InformationSet<TemplatedState> InformationSet;

	enum class Move : signed char { rock, paper, scissors, spock, lizard };
	typedef std::array<double, 2> Reward;

	static int num_players() { return 2; }

	static std::string move_name(Move move) {
		if (move == Move::rock) {
			return "rock    ";
		} else if (move == Move::paper) {
			return "paper   ";
		} else if (move == Move::scissors) {
			return "scissors";
		} else if (move == Move::spock) {
			return "Spock";
		} else if (move == Move::lizard) {
			return "lizard";
		} else {
			throw std::invalid_argument("Invalid move.");
		}
	}

	// Returns a vector of all possible initial states (with equal probability).
	// This is needed for solving the game exactly.
	static std::vector<TemplatedState> all_initial_states() { return {TemplatedState()}; }

	bool terminal() const { return game_round >= 2; }

	int player() const { return game_round % 2; }

	Reward reward() const {
		minimum_core_assert(terminal());
		return {double(game_result), double(-game_result)};
	}

	std::vector<Move> possible_moves() const {
		minimum_core_assert(!terminal());

		if (include_spock) {
			if (player() == 0) {
				// Remove the Spock option for assymmetry.
				return {Move::rock, Move::paper, Move::scissors, /*Move::spock,*/ Move::lizard};
			} else {
				return {Move::rock, Move::paper, Move::scissors, Move::spock, Move::lizard};
			}
		} else {
			return {Move::rock, Move::paper, Move::scissors};
		}
	}

	void move(Move move) {
		minimum_core_assert(!terminal());
		if (game_round == 0) {
			player0_move = move;
		} else {
			player1_move = move;
			if (check_win(player0_move, player1_move)) {
				game_result = 1;
			} else if (check_win(player1_move, player0_move)) {
				game_result = -1;
			} else {
				game_result = 0;
			}
		}
		game_round++;
	}

	bool check_win(Move move1, Move move2) {
		if (move1 == Move::rock) {
			return move2 == Move::scissors || move2 == Move::lizard;
		} else if (move1 == Move::paper) {
			return move2 == Move::rock || move2 == Move::spock;
		} else if (move1 == Move::scissors) {
			return move2 == Move::paper || move2 == Move::lizard;
		} else if (move1 == Move::spock) {
			return move2 == Move::rock || move2 == Move::scissors;
		} else if (move1 == Move::lizard) {
			return move2 == Move::paper || move2 == Move::spock;
		}
		return false;
	}

   private:
	Move player0_move = Move::rock;
	Move player1_move = Move::rock;

	signed char game_round = 0;
	signed char game_result = 0;
};

using State = TemplatedState<false>;
using SpockState = TemplatedState<true>;

// The players know nothing in Rock, Paper, Scissors.
//
// But we still need to differentiate between the null information
// set created by the default constructor.
template <typename State>
class InformationSet {
   public:
	typedef typename State::Move Move;
	InformationSet() : is_null(true) {}
	InformationSet(const State& state) : is_null(false) {}
	bool operator==(const InformationSet& other) const { return is_null == other.is_null; }
	std::string str() const { return "<No information>"; }

	bool is_null = false;
};

}  // namespace rock_paper_scissors

namespace std {
template <typename State>
struct hash<rock_paper_scissors::InformationSet<State>> {
	size_t operator()(rock_paper_scissors::InformationSet<State> const& is) const { return 0; }
};

template <>
struct hash<typename rock_paper_scissors::State::Move> {
	std::hash<signed char> hasher;
	size_t operator()(typename rock_paper_scissors::State::Move const& move) const {
		return hasher((unsigned char)move);
	}
};

template <>
struct hash<typename rock_paper_scissors::SpockState::Move> {
	std::hash<signed char> hasher;
	size_t operator()(typename rock_paper_scissors::SpockState::Move const& move) const {
		return hasher((unsigned char)move);
	}
};
}  // namespace std

#endif
