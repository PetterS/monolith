#pragma once
// Algorithm from
//
// [1] Fast Algorithms for Finding Randomized Strategies in Game Trees
// http://www.maths.lse.ac.uk/personal/stengel/TEXTE/stoc.pdf

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
#include <minimum/core/range.h>
#include <minimum/core/string.h>
#include <minimum/linear/ip.h>
#include <minimum/linear/solver.h>
using minimum::core::range;
using minimum::core::Timer;
using minimum::core::to_string;

namespace minimum {
namespace ai {
namespace equilibrium_sequential {

template <typename State>
struct Action {
	typename State::InformationSet from;
	typename State::Move move;
	bool empty;

	Action() : empty(true) {}
	Action(typename State::InformationSet from_, typename State::Move move_)
	    : from(from_), move(move_), empty(false) {
		typename State::InformationSet null;
		minimum::core::check(!(from == null),
		                     "The \"null\" information set can not be equal to any other!");
	}

	bool operator==(const Action<State>& rhs) const {
		if (empty != rhs.empty) {
			return false;
		} else if (empty) {
			return true;
		}
		return from == rhs.from && move == rhs.move;
	}

	struct Hash {
		typedef Action<State> argument_type;
		typedef size_t result_type;
		std::hash<typename State::InformationSet> hash1;
		std::hash<typename State::Move> hash2;
		result_type operator()(const argument_type& arg) const {
			if (arg.empty) {
				return 0xdead;
			}
			return hash1(arg.from) ^ hash2(arg.move);
		}
	};

	struct PairHash {
		typedef std::pair<Action<State>, Action<State>> argument_type;
		typedef size_t result_type;
		Hash hash;
		result_type operator()(const argument_type& arg) const {
			return hash(arg.first) ^ (hash(arg.second) >> 1);
		}
	};
};

template <typename State>
std::ostream& operator<<(std::ostream& out, const Action<State>& action) {
	if (action.empty) {
		out << 0;
	} else {
		out << State::move_name(action.move) << " when " << action.from.str();
		// out << State::move_name(action.move);
	}
	return out;
}

template <typename State>
struct RecurseInfo {
	std::vector<int> parent[2];
	std::unordered_map<Action<State>, int, typename Action<State>::Hash> player_sequences[2];

	std::unordered_map<std::pair<Action<State>, Action<State>>,
	                   double,
	                   typename Action<State>::PairHash>
	    game_values;

	double probability_of_each_initial_state = 0;

	// Used to detect games of imperfect recall, which do not work.
	std::unordered_map<typename State::InformationSet, typename State::InformationSet>
	    loop_detection;
};

template <typename State, typename InformationSet>
void recurse_tree(State state,
                  std::vector<std::vector<Action<State>>>* current_sets,
                  RecurseInfo<State>* recurse_info) {
	using namespace std;

	if (state.terminal()) {
		// cout << "Terminal.\n";
		// cout << "-- P0: " << to_string(current_sets->at(0)) << "\n";
		// cout << "-- P1: " << to_string(current_sets->at(1)) << "\n";
		for (int i = 0; i < 2; ++i) {
			auto& sequences = recurse_info->player_sequences[i];
			int current_parent = -1;
			for (const auto& action : current_sets->at(i)) {
				auto itr = sequences.find(action);
				if (itr == sequences.end()) {
					int new_index = int(recurse_info->parent[i].size());
					sequences[action] = new_index;
					recurse_info->parent[i].emplace_back(current_parent);
					current_parent = new_index;
				} else {
					current_parent = itr->second;
				}
			}
		}
		auto value = state.reward()[0] * recurse_info->probability_of_each_initial_state;
		// cout << "-- Reward: " << value << "\n";
		recurse_info
		    ->game_values[make_pair(current_sets->at(0).back(), current_sets->at(1).back())] +=
		    value;
		return;
	}

	int i = state.player();
	minimum_core_assert(i == 0 || i == 1);

	auto moves = state.possible_moves();
	InformationSet is(state);

	// Check that this information state has a unique predecessor (for this player).
	// Otherwise, we have a game of imperfect recall, resulting in bugs that can
	// be hard to find.
	const auto& prev_is = current_sets->at(i).back().from;
	if (recurse_info->loop_detection.count(is) > 0) {
		minimum::core::check(recurse_info->loop_detection[is] == prev_is,
		                     "The information sets for player ",
		                     i,
		                     " do not form a tree.\n"
		                     "This seems to be a game of imperfect recall.\n\n"
		                     "Current information set: ",
		                     is.str(),
		                     "\n"
		                     "Predecessor 1: ",
		                     prev_is.str(),
		                     "\n"
		                     "Predecessor 2: ",
		                     recurse_info->loop_detection[is].str(),
		                     "\n");
	}
	recurse_info->loop_detection.emplace(is, prev_is);

	for (auto move : moves) {
		State new_state(state);
		new_state.move(move);
		current_sets->at(i).emplace_back(is, move);
		recurse_tree<State, InformationSet>(new_state, current_sets, recurse_info);
		current_sets->at(i).pop_back();
	}
}

template <typename State, typename InformationSet>
void explore_state(State state, RecurseInfo<State>* recurse_info) {
	using namespace std;
	vector<vector<Action<State>>> current_sets(State::num_players());
	for (auto& actions : current_sets) {
		actions.emplace_back();
	}
	recurse_tree<State, InformationSet>(state, &current_sets, recurse_info);
}

template <typename State>
void all_sequences(RecurseInfo<State>* recurse_info) {
	using namespace std;
	typedef typename State::InformationSet InformationSet;

	// All possible (hidden) states the game can start in. They occur
	// with equal probability.
	auto all_initial_states = State::all_initial_states();
	recurse_info->probability_of_each_initial_state = 1.0 / double(all_initial_states.size());

	// Assume a two-player game.
	minimum_core_assert(State::num_players() == 2);

	// Compute all sequences each player has.
	for (const auto& state : all_initial_states) {
		explore_state<State, InformationSet>(state, recurse_info);
	}
}

template <typename State>
void compute_matrices(const RecurseInfo<State> recurse_info,
                      Eigen::SparseMatrix<double>* E,
                      std::vector<double>* e,
                      Eigen::SparseMatrix<double>* F,
                      std::vector<double>* f) {
	using namespace std;
	typedef typename State::InformationSet InformationSet;

	// Assign all information sets an index for use in the matrices.
	unordered_map<InformationSet, int> player_sets[2];
	for (int i = 0; i < 2; ++i) {
		player_sets[i][InformationSet{}] = 0;

		for (const auto& elem : recurse_info.player_sequences[i]) {
			if (!elem.first.empty) {
				auto& from_set = elem.first.from;
				if (player_sets[i].count(from_set) == 0) {
					int new_index = int(player_sets[i].size());
					player_sets[i][from_set] = new_index;
				}
			}
		}
	}

	vector<Eigen::Triplet<double>> Es[2];

	for (int i = 0; i < 2; ++i) {
		// Whether we have added the -1 coefficient for the parent for the particular row.
		vector<bool> have_parent(player_sets[i].size(), false);

		int empty_row = player_sets[i][InformationSet{}];
		int empty_col = recurse_info.player_sequences[i].at(Action<State>{});
		Es[i].emplace_back(empty_row, empty_col, 1.0);

		for (auto& elem : recurse_info.player_sequences[i]) {
			auto& action = elem.first;
			int action_col = elem.second;

			if (!action.empty) {
				int from_information_set = player_sets[i].at(action.from);
				minimum_core_assert(from_information_set != empty_row);
				Es[i].emplace_back(from_information_set, action_col, 1.0);

				int parent = recurse_info.parent[i].at(action_col);
				if (parent >= 0 && !have_parent.at(from_information_set)) {
					Es[i].emplace_back(from_information_set, parent, -1.0);
					have_parent.at(from_information_set) = true;
				}
			}
		}
	}

	E->resize(int(player_sets[0].size()), int(recurse_info.player_sequences[0].size()));
	E->setFromTriplets(begin(Es[0]), end(Es[0]));
	e->clear();
	e->resize(player_sets[0].size(), 0.0);
	e->at(player_sets[0][InformationSet{}]) = 1;

	F->resize(int(player_sets[1].size()), int(recurse_info.player_sequences[1].size()));
	F->setFromTriplets(begin(Es[1]), end(Es[1]));
	f->clear();
	f->resize(player_sets[1].size(), 0.0);
	f->at(player_sets[1][InformationSet{}]) = 1;
}

template <typename State>
void compute_matrix_A(const RecurseInfo<State> recurse_info, Eigen::SparseMatrix<double>* A) {
	std::vector<Eigen::Triplet<double>> As;
	for (auto& elem : recurse_info.game_values) {
		auto& action_pair = elem.first;
		double value = elem.second;

		int row = recurse_info.player_sequences[0].at(action_pair.first);
		int col = recurse_info.player_sequences[1].at(action_pair.second);
		As.emplace_back(row, col, value);
	}
	A->resize(int(recurse_info.player_sequences[0].size()),
	          int(recurse_info.player_sequences[1].size()));
	A->setFromTriplets(begin(As), end(As));
}

template <typename State>
struct Result {
	double value = 0;

	std::unordered_map<typename State::InformationSet,
	                   std::vector<std::pair<typename State::Move, double>>>
	    player_strategies[2];
};

template <typename State>
void extract_solution(const RecurseInfo<State>& recurse_info,
                      const Eigen::VectorXd& solution,
                      int player,
                      Result<State>* result) {
	for (auto& elem : recurse_info.player_sequences[player]) {
		auto& action = elem.first;
		if (!action.empty) {
			auto index = elem.second;
			double probability = solution(index);

			int parent = recurse_info.parent[player].at(index);
			minimum_core_assert(parent >= 0);
			double parent_probability = solution(parent);
			probability /= std::max(parent_probability, 1e-9);

			auto& information_set = action.from;
			auto& move = action.move;
			result->player_strategies[player][information_set].emplace_back(move, probability);
		}
	}
}

template <typename State>
void print_result(const Result<State>& result) {
	using namespace std;
	cout << "Game value is " << result.value << " for player 0.\n";
	cout << "Player 0 has " << result.player_strategies[0].size() << " information states.\n";
	cout << "Player 1 has " << result.player_strategies[1].size() << " information states.\n";
	cout << "\n";

	for (int i = 0; i < 2; ++i) {
		cout << "Player " << i << "\n";
		for (auto& elem : result.player_strategies[i]) {
			auto& information_set = elem.first;
			auto& moves = elem.second;
			if (moves.empty() || moves.back().second != moves.back().second) {
				continue;
			}
			cout << "-- At information set " << information_set.str() << ":\n";
			for (auto& move : moves) {
				cout << "---- Do " << State::move_name(move.first) << " with probability "
				     << move.second << "\n";
			}
		}
		cout << "\n";
	}
}

namespace {
void prefer_integer_solution(minimum::linear::IP* ip,
                             const minimum::linear::Sum& objective,
                             const std::vector<minimum::linear::Variable>& x) {
	// Find a solution as little fractional as possible.
	if (x.size() < 100) {
		ip->clear_objective();
		ip->add_constraint(objective == objective.value());
		for (auto i : range(x.size())) {
			auto fractional = ip->add_variable(minimum::linear::IP::Real);
			ip->add_bounds(0, fractional, 1);
			auto boolean = ip->add_boolean();
			ip->add_constraint(x[i] == fractional + boolean);
			auto fractional_used = ip->add_boolean();
			ip->add_constraint(fractional <= fractional_used);
			ip->add_objective(fractional_used);
		}
		minimum::core::check(solve(ip), "Could not solve MIP.");
	}
}
}  // namespace

template <typename State>
Result<State> compute() {
	using namespace std;
	using namespace minimum::linear;

	RecurseInfo<State> recurse_info;
	all_sequences(&recurse_info);

	// for (auto i : {0, 1}) {
	// 	clog << "Player " << i << " has " << recurse_info.player_sequences[i].size() << "
	// actions.\n"; 	clog << "-- " << to_string(recurse_info.player_sequences[i]) << "\n";
	// }

	Eigen::SparseMatrix<double> E;
	vector<double> e;
	Eigen::SparseMatrix<double> F;
	vector<double> f;
	compute_matrices(recurse_info, &E, &e, &F, &f);

	Eigen::SparseMatrix<double> A;
	compute_matrix_A(recurse_info, &A);

	Eigen::VectorXd x_sol(A.rows());
	Eigen::VectorXd y_sol(A.cols());

	// clog << E.toDense() << "\n\n";
	// clog << F.toDense() << "\n\n";
	// clog << A.toDense() << "\n\n";

	Result<State> result;
	{
		Timer t("Solving LP.");
		IP ip;
		auto y = ip.add_vector(A.cols(), IP::Real);
		auto p = ip.add_vector(E.rows(), IP::Real);

		Sum eTp = 0;
		for (int i = 0; i < E.rows(); ++i) {
			eTp += e[i] * p[i];
		}
		ip.add_objective(eTp);

		// -Ay + E^Tp ≥ 0.
		vector<DualVariable> x;
		for (int i = 0; i < A.rows(); ++i) {
			Sum row = 0;
			for (int j = 0; j < A.cols(); ++j) {
				row += -A.coeff(i, j) * y[j];
			}
			for (int j = 0; j < E.rows(); ++j) {
				row += E.coeff(j, i) * p[j];
			}
			x.push_back(ip.add_constraint(row >= 0));
		}

		// -Fy = -f.
		for (int i = 0; i < F.rows(); ++i) {
			Sum row = 0;
			for (int j = 0; j < F.cols(); ++j) {
				row += -F.coeff(i, j) * y[j];
			}
			ip.add_constraint(row == -f[i]);
		}

		// y ≥ 0.
		for (int j = 0; j < A.cols(); ++j) {
			ip.add_constraint(y[j] >= 0);
		}

		minimum::core::check(solve(&ip), "Could not solve linear program.");
		result.value = eTp.value();

		prefer_integer_solution(&ip, eTp, y);
		for (int j = 0; j < A.cols(); ++j) {
			y_sol(j) = y.at(j).value();
		}
		extract_solution(recurse_info, y_sol, 1, &result);
		t.OK();
	}

	// First verification: Problem (9) from [1].
	{
		Timer t("Verification (9).");
		IP ip;
		auto x = ip.add_vector(A.rows(), IP::Real);
		auto q = ip.add_vector(F.rows(), IP::Real);

		Sum qTf = 0;
		for (int i = 0; i < F.rows(); ++i) {
			qTf += f[i] * q[i];
		}
		ip.add_objective(qTf);

		// x^T(-A) - q^TF ≤ 0.
		for (int j = 0; j < A.cols(); ++j) {
			Sum col = 0;
			for (int i = 0; i < A.rows(); ++i) {
				col += -A.coeff(i, j) * x[i];
			}
			for (int i = 0; i < F.rows(); ++i) {
				col += -F.coeff(i, j) * q[i];
			}
			ip.add_constraint(col <= 0);
		}

		// Ex =	e.
		for (int i = 0; i < E.rows(); ++i) {
			Sum row = 0;
			for (int j = 0; j < E.cols(); ++j) {
				row += E.coeff(i, j) * x[j];
			}
			ip.add_constraint(row == e[i]);
		}

		// x ≥ 0.
		for (int j = 0; j < A.rows(); ++j) {
			ip.add_constraint(x[j] >= 0);
		}

		minimum::core::check(solve(&ip), "Could not solve linear program (9).");
		minimum_core_assert(abs(qTf.value() + result.value) <= 1e-6);

		for (int i = 0; i < A.rows(); ++i) {
			x_sol(i) = x.at(i).value();
		}
		prefer_integer_solution(&ip, qTf, x);
		extract_solution(recurse_info, x_sol, 0, &result);

		t.OK();
	}

	// Second verification. Best response for player 0 to player 1’s strategy.
	// Problem (5) from [1].
	{
		Timer t("Verification (5).");
		IP ip;

		auto x = ip.add_vector(A.rows(), IP::Real);
		Eigen::VectorXd Ay = A * y_sol;

		Sum xAy = 0;
		for (int i = 0; i < A.rows(); ++i) {
			xAy += Ay(i) * x[i];
			ip.add_constraint(x[i] >= 0);
		}
		ip.add_objective(-xAy);

		// Ex =	e.
		for (int i = 0; i < E.rows(); ++i) {
			Sum row = 0;
			for (int j = 0; j < E.cols(); ++j) {
				row += E.coeff(i, j) * x[j];
			}
			ip.add_constraint(row == e[i]);
		}

		minimum::core::check(solve(&ip), "Could not solve linear program (5).");
		minimum_core_assert(abs(xAy.value() - result.value) <= 1e-6);
		t.OK();
	}

	// Third verification. Best response for player 1 to player 0’s strategy.
	// Problem (2) from [1].
	{
		Timer t("Verification (2).");
		IP ip;

		auto y = ip.add_vector(A.cols(), IP::Real);
		Eigen::VectorXd ATx = A.transpose() * x_sol;

		Sum yATx = 0;
		for (int j = 0; j < A.cols(); ++j) {
			yATx += ATx(j) * y[j];
			ip.add_constraint(y[j] >= 0);
		}
		ip.add_objective(yATx);

		// Fy = f.
		for (int i = 0; i < F.rows(); ++i) {
			Sum row = 0;
			for (int j = 0; j < F.cols(); ++j) {
				row += F.coeff(i, j) * y[j];
			}
			ip.add_constraint(row == f[i]);
		}

		minimum::core::check(solve(&ip), "Could not solve linear program (2).");
		minimum_core_assert(abs(yATx.value() - result.value) <= 1e-6);
		t.OK();
	}

	return result;
}
}  // namespace equilibrium_sequential
}  // namespace ai
}  // namespace minimum
