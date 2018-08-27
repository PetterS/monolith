#pragma once

// Computes a mixed Nash equilirium for a game with imperfect information
// using linear programming.
//
// A game with imperfect information is defined by two classes:
//
//  • class State – containing all information about the current game state.
//  • ⁠class InformationSet – containing all information that a player knows.
//                           A player can not distingusih between equal
//                           information states.
//
// See kuhn_poker.h or rock_paper_scissors.h for simple examples of how
// these classes define a game. The key operation is the constructor
//
//    InformationSet::InformationSet(const State& state);
//
// which defines how much of the game state players actually observe.

#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <minimum/core/check.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using minimum::core::to_string;

namespace minimum {
namespace ai {
namespace equilibrium {

template <typename InformationSet>
struct Result {
	typedef typename InformationSet::Move Move;
	typedef std::vector<Move> Moves;
	// A strategy is a mapping from every information state to a move.
	typedef std::unordered_map<InformationSet, Move> Strategy;

	// Value of the game for the first player.
	double value = 0;

	// player_strategies[0] contains all strategies with probabilities for
	// player 0 etc.
	std::vector<std::vector<std::pair<double, Strategy>>> player_strategies;
};

// Computes the game value given pure strategies for both players.
template <typename State, typename Strategy>
double game_value(State state, const Strategy& strategy0, const Strategy& strategy1) {
	while (!state.terminal()) {
		typename State::InformationSet set(state);
		typename State::Move move;
		if (state.player() == 0) {
			auto move_itr = strategy0.find(set);
			minimum_core_assert(move_itr != end(strategy0));
			move = move_itr->second;
		} else {
			auto move_itr = strategy1.find(set);
			minimum_core_assert(move_itr != end(strategy1));
			move = move_itr->second;
		}
		state.move(move);
	}
	auto reward = state.reward();
	minimum_core_assert(reward[0] == -reward[1]);  // Zero-sum game.
	return reward[0];
}

template <typename State, typename InformationSet, typename Strategy>
void create_game_matrices(const std::vector<std::vector<Strategy>> strategies,
                          std::vector<std::vector<double>>* M,
                          std::vector<std::vector<double>>* MT,
                          double* value_offset) {
	auto all_initial_states = State::all_initial_states();

	*value_offset = 0;

	for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
		M->emplace_back();
		for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
			M->at(s1).emplace_back(0);

			// For these two strategies, go through all possible starting
			// game states.
			double value = 0;
			double games = 0;
			for (const auto& state : all_initial_states) {
				value += game_value(state, strategies[0][s0], strategies[1][s1]);
				games += 1;
			}
			value /= games;
			(*M)[s1][s0] = value;
			*value_offset = std::max(*value_offset, value);
		}
	}

	for (auto& row : *M) {
		for (auto& value : row) {
			value += *value_offset;
		}
	}

	for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
		(*MT).emplace_back();
		for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
			(*MT)[s0].emplace_back(0);
			(*MT)[s0][s1] = -(*M)[s1][s0] + 2 * (*value_offset);
		}
	}
}

template <typename State, typename InformationSet, typename Strategy, typename Decisions>
std::vector<std::vector<Strategy>> enumerate_strategies(int num_players,
                                                        const Decisions& decisions) {
	using namespace std;

	// Each player has a vector of strategies she can choose from.
	vector<vector<Strategy>> strategies(num_players);

	// Compute all possible strategies
	for (int player = 0; player < 2; ++player) {
		unordered_map<InformationSet, size_t> current_pos;
		while (true) {
			Strategy strategy;
			for (auto& entry : decisions[player]) {
				auto& set = entry.first;
				auto& moves = entry.second;
				auto& index = current_pos[set];
				strategy[set] = moves[index];
			}
			strategies[player].emplace_back(move(strategy));

			auto itr = decisions[player].begin();
			current_pos[itr->first]++;
			bool done = false;
			while (current_pos[itr->first] >= itr->second.size()) {
				current_pos[itr->first] = 0;
				itr++;
				if (itr == decisions[player].end()) {
					done = true;
					break;
				}
				current_pos[itr->first]++;
			}
			if (done) {
				break;
			}
		}
		cout << "Player " << player << " has " << strategies[player].size() << " strategies.\n";
	}
	return strategies;
}

template <typename State, typename InformationSet, typename Sets>
void recurse_tree(State state,
                  std::vector<typename Result<InformationSet>::Moves>* current_moves,
                  Sets* sets) {
	using namespace std;

	if (state.terminal()) {
		return;
	}

	InformationSet set(state);
	auto moves = state.possible_moves();
	int i = state.player();
	for (auto move : moves) {
		State new_state(state);
		new_state.move(move);
		current_moves->at(i).push_back(move);
		recurse_tree<State, InformationSet>(new_state, current_moves, sets);
		current_moves->at(i).pop_back();
	}

	sets->at(i).insert(make_pair(set, moves));
}

template <typename State, typename InformationSet, typename Sets>
void explore_state(State state, Sets* sets) {
	using namespace std;
	vector<typename Result<InformationSet>::Moves> current_moves(State::num_players());
	recurse_tree<State, InformationSet>(state, &current_moves, sets);
}

template <typename State>
std::vector<std::unordered_map<typename State::InformationSet, std::vector<typename State::Move>>>
all_decisions() {
	using namespace std;
	typedef typename State::InformationSet InformationSet;

	// All possible (hidden) states the game can start in. They occur
	// with equal probability.
	auto all_initial_states = State::all_initial_states();

	// Assume a two-player game.
	minimum_core_assert(State::num_players() == 2);

	// Compute all decisions each player has.
	vector<unordered_map<InformationSet, vector<typename State::Move>>> decisions(
	    State::num_players());
	for (const auto& state : all_initial_states) {
		explore_state<State, InformationSet>(state, &decisions);
	}

	return decisions;
}

template <typename State, typename Matrix, typename Strategies>
void print_game_matrix(const Matrix& M, double value_offset, const Strategies& strategies) {
	using namespace std;

	if (max(strategies[0].size(), strategies[1].size()) > 10) {
		return;
	}

	cout << "       ";
	for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
		cout << "(";
		for (auto&& elem : strategies[1][s1]) {
			cout << State::move_name(elem.second);
		}
		cout << ") ";
	}
	cout << "\n";
	for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
		cout << "(";
		for (auto&& elem : strategies[0][s0]) {
			cout << State::move_name(elem.second);
		}
		cout << ")  ";
		for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
			cout << setw(6) << right << M.at(s1).at(s0) - value_offset;
		}
		cout << "\n";
	}
}

template <typename State>
Result<typename State::InformationSet> compute() {
	using namespace std;
	using namespace minimum::linear;
	typedef typename State::InformationSet InformationSet;

	// All possible (hidden) states the game can start in. They occur
	// with equal probability.
	auto all_initial_states = State::all_initial_states();

	// Compute all decisions each player has.
	auto decisions = all_decisions<State>();

	// Player 0 always folds with the best card in the end.
	// decisions[0][InformationSet(State(3, 1).move(State::Move::check).move(State::Move::bet))] =
	// {State::Move::fold};

	// Player 0 always checks in the beginning.
	// decisions[0][InformationSet(State(3, 1))] = {State::Move::check};
	// decisions[0][InformationSet(State(2, 1))] = {State::Move::check};
	// decisions[0][InformationSet(State(1, 2))] = {State::Move::check};

	// Print all decisions.
	for (int player = 0; player < 2; ++player) {
		cout << "Player " << player << '\n';
		for (auto& entry : decisions[player]) {
			auto& set = entry.first;
			auto& moves = entry.second;

			vector<string> string_moves;
			for (auto&& move : moves) {
				string_moves.emplace_back(State::move_name(move));
			}
			cout << "Decision: " << set.str() << " with moves " << to_string(string_moves) << "\n";
		}
		cout << '\n';
	}

	typedef typename Result<InformationSet>::Strategy Strategy;
	auto strategies =
	    enumerate_strategies<State, InformationSet, Strategy>(State::num_players(), decisions);

	// Create game matrices (notation from http://en.wikipedia.org/wiki/Zero-sum_game),
	// containing the payoffs for every combination of strategies.
	vector<vector<double>> M, MT;
	double value_offset = 0;
	create_game_matrices<State, InformationSet>(strategies, &M, &MT, &value_offset);

	print_game_matrix<State>(M, value_offset, strategies);

	Result<InformationSet> result;
	result.player_strategies.resize(2);

	// Compute the optimal strategy for the first player.
	{
		IP ip;
		auto u = ip.add_vector(strategies[0].size(), IP::Real);
		for (auto& ui : u) {
			ip.add_constraint(ui >= 0);
			ip.add_objective(ui);
		}
		for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
			Sum sum = 0;
			for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
				sum += M[s1][s0] * u[s0];
			}
			ip.add_constraint(sum >= 1);
		}
		minimum_core_assert(solve(&ip));

		double usum = 0;
		for (auto& ui : u) {
			usum += ui.value();
		}
		double value = 1.0 / usum;
		result.value = value - value_offset;

		for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
			if (u[s0].value() > 1e-9) {
				double p = u[s0].value() / usum;
				result.player_strategies[0].emplace_back(p, strategies[0][s0]);
			}
		}
	}

	// Compute the optimal strategy for the second player (because we don’t have the duals).
	{
		IP ip;
		auto u = ip.add_vector(strategies[1].size(), IP::Real);
		for (auto& ui : u) {
			ip.add_constraint(ui >= 0);
			ip.add_objective(ui);
		}
		for (int s0 = 0; s0 < strategies[0].size(); ++s0) {
			Sum sum = 0;
			for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
				sum += MT[s0][s1] * u[s1];
			}
			ip.add_constraint(sum >= 1);
		}
		minimum_core_assert(solve(&ip));

		double usum = 0;
		for (auto& ui : u) {
			usum += ui.value();
		}
		double value = 1.0 / usum;
		minimum_core_assert(result.value + (value - value_offset) <= 1e-6);

		for (int s1 = 0; s1 < strategies[1].size(); ++s1) {
			if (u[s1].value() > 1e-9) {
				double p = u[s1].value() / usum;
				result.player_strategies[1].emplace_back(p, strategies[1][s1]);
			}
		}
	}

	// Verify that no pure strategy for player 0 gives a better result.
	for (const auto& strategy : strategies[0]) {
		double value = 0;
		for (const auto& strategy_pair1 : result.player_strategies[1]) {
			for (const auto& state : all_initial_states) {
				auto p = strategy_pair1.first;
				value += p * game_value(state, strategy, strategy_pair1.second);
			}
		}
		value /= all_initial_states.size();
		minimum_core_assert(value <= result.value + 1e-6);
	}

	// Verify the equivalent for player 1.
	for (const auto& strategy : strategies[1]) {
		double value = 0;
		for (const auto& strategy_pair0 : result.player_strategies[0]) {
			for (const auto& state : all_initial_states) {
				auto p = strategy_pair0.first;
				value += p * game_value(state, strategy_pair0.second, strategy);
			}
		}
		value /= all_initial_states.size();
		minimum_core_assert(value >= result.value - 1e-6);
	}

	return result;
}
}  // namespace equilibrium
}  // namespace ai
}  // namespace minimum
