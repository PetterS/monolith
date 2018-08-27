#pragma once
// Game from
// Fast Algorithms for Finding Randomized Strategies in Game Trees
// http://www.maths.lse.ac.uk/personal/stengel/TEXTE/stoc.pdf

#include <array>
#include <functional>
#include <random>
#include <sstream>

#include <minimum/core/check.h>

namespace stengel {

class InformationSet;

class State {
   public:
	friend class InformationSet;
	typedef stengel::InformationSet InformationSet;

	enum class Move : signed char { l, r, L, R, c, d, p, q, s, t };
	typedef std::array<double, 2> Reward;

	static int num_players() { return 2; }

	static std::string move_name(Move move) {
		if (move == Move::l) {
			return "l";
		} else if (move == Move::r) {
			return "r";
		} else if (move == Move::L) {
			return "L";
		} else if (move == Move::R) {
			return "R";
		} else if (move == Move::c) {
			return "c";
		} else if (move == Move::d) {
			return "d";
		} else if (move == Move::p) {
			return "p";
		} else if (move == Move::q) {
			return "q";
		} else if (move == Move::s) {
			return "s";
		} else if (move == Move::t) {
			return "t";
		} else {
			throw std::invalid_argument("Invalid move.");
		}
	}

	// Returns a vector of all possible initial states (with equal probability).
	// This is needed for solving the game exactly.
	static std::vector<State> all_initial_states() {
		return {State(1), State(2), State(3), State(4), State(4)};
	}

	State(int initial_position) {
		minimum_core_assert(1 <= initial_position && initial_position <= 4);
		this->state = initial_position;
		if (initial_position == 4) {
			game_result = 5;
			state = 100;
		}
	}

	bool terminal() const { return state >= 9; }

	int player() const {
		if (4 <= state && state <= 6) {
			return 1;
		} else {
			return 0;
		}
	}

	Reward reward() const {
		minimum_core_assert(terminal());
		return {game_result, -game_result};
	}

	std::vector<Move> possible_moves() const {
		minimum_core_assert(!terminal());
		minimum_core_assert(1 <= state && state <= 8);
		if (state == 1) {
			return {Move::l, Move::r};
		} else if (2 <= state && state <= 3) {
			return {Move::c, Move::d};
		} else if (state == 4) {
			return {Move::p, Move::q};
		} else if (5 <= state && state <= 6) {
			return {Move::s, Move::t};
		} else {
			return {Move::L, Move::R};
		}
	}

	void move(Move move) {
		minimum_core_assert(!terminal());
		minimum_core_assert(1 <= state && state <= 8);
		if (state == 1) {
			minimum_core_assert(move == Move::l || move == Move::r);
			if (move == Move::l) {
				game_result = 5;
				state = 100;
			} else {
				state = 4;
			}
		} else if (state == 2) {
			minimum_core_assert(move == Move::c || move == Move::d);
			if (move == Move::c) {
				game_result = 10;
			} else {
				game_result = 20;
			}
			state = 100;
		} else if (state == 3) {
			minimum_core_assert(move == Move::c || move == Move::d);
			if (move == Move::c) {
				state = 5;
			} else {
				state = 6;
			}
		} else if (state == 4) {
			minimum_core_assert(move == Move::p || move == Move::q);
			if (move == Move::p) {
				state = 7;
			} else {
				state = 8;
			}
		} else if (state == 5) {
			minimum_core_assert(move == Move::s || move == Move::t);
			if (move == Move::s) {
				game_result = 20;
			} else {
				game_result = 50;
			}
			state = 100;
		} else if (state == 6) {
			minimum_core_assert(move == Move::s || move == Move::t);
			if (move == Move::s) {
				game_result = 30;
			} else {
				game_result = 15;
			}
			state = 100;
		} else if (state == 7) {
			minimum_core_assert(move == Move::L || move == Move::R);
			if (move == Move::L) {
				game_result = 10;
			} else {
				game_result = 15;
			}
			state = 100;
		} else if (state == 8) {
			minimum_core_assert(move == Move::L || move == Move::R);
			if (move == Move::L) {
				game_result = 20;
			} else {
				game_result = -5;
			}
			state = 100;
		}
	}

   private:
	int state = 0;
	double game_result = -1e9;
};

class InformationSet {
   public:
	typedef State::Move Move;
	InformationSet() {}
	InformationSet(const State& state) {
		if (state.terminal()) {
			information_set = 100;
			return;
		}

		auto id = state.state;
		minimum_core_assert(1 <= id && id <= 8);
		if (id == 1) {
			information_set = 1;
		} else if (2 <= id && id <= 3) {
			information_set = 2;
		} else if (id == 4) {
			information_set = 3;
		} else if (5 <= id && id <= 6) {
			information_set = 4;
		} else if (7 <= id && id <= 8) {
			information_set = 5;
		}
	}
	bool operator==(const InformationSet& other) const {
		return information_set == other.information_set;
	}
	std::string str() const { return "IS-" + std::to_string(information_set); }

	int hash() const { return information_set; }

   private:
	int information_set = -1;
};

}  // namespace stengel

namespace std {
template <>
struct hash<stengel::InformationSet> {
	size_t operator()(const stengel::InformationSet& is) const { return is.hash(); }
};

template <>
struct hash<typename stengel::State::Move> {
	std::hash<signed char> hasher;
	size_t operator()(const stengel::State::Move& move) const {
		return hasher((unsigned char)move);
	}
};
}  // namespace std
