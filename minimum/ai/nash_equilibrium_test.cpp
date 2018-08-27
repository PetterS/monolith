#include <catch.hpp>

#include <minimum/ai/compute_equilibrium.h>
#include <minimum/ai/compute_equilibrium_sequential.h>
#include <minimum/ai/imperfect_games/kuhn_poker.h>
#include <minimum/ai/imperfect_games/leduc_holdem.h>
#include <minimum/ai/imperfect_games/rock_paper_scissors.h>
#include <minimum/ai/imperfect_games/stengel.h>

using namespace std;
using namespace minimum::ai;

TEST_CASE("rock_paper_scissors") {
	auto result = equilibrium::compute<rock_paper_scissors::State>();
	REQUIRE(result.player_strategies.size() == 2);
	REQUIRE(result.player_strategies[0].size() == 3);
	REQUIRE(result.player_strategies[1].size() == 3);
	for (int player : {0, 1}) {
		for (auto pure : result.player_strategies[0]) {
			CHECK(Approx(pure.first) == 1.0 / 3.0);
		}
	}
	CHECK(std::abs(result.value) <= 1e-6);
}

TEST_CASE("kuhn_poker") {
	typedef kuhn_poker::State<> State;
	typedef kuhn_poker::InformationSet<State> IS;

	auto result = equilibrium::compute<State>();
	REQUIRE(result.player_strategies.size() == 2);
	CHECK(Approx(result.value) == -1.0 / 18.0);

	auto result2 = equilibrium_sequential::compute<State>();
	CHECK(Approx(result2.value) == -1.0 / 18.0);

	// Second player has a unique strategy.
	// See “Effective Short-Term Opponent Exploitation in Simplified Poker.”
	{
		State state(2, 1);
		state.move(State::Move::check);
		auto moves = result2.player_strategies[1].at(IS(state));
		bool found = false;
		for (auto move : moves) {
			if (move.first == State::Move::bet) {
				CHECK(Approx(move.second) == 1.0 / 3.0);
				found = true;
			}
		}
		CHECK(found);
	}
	{
		State state(3, 2);
		state.move(State::Move::bet);
		auto moves = result2.player_strategies[1].at(IS(state));
		bool found = false;
		for (auto move : moves) {
			if (move.first == State::Move::call) {
				CHECK(Approx(move.second) == 1.0 / 3.0);
				found = true;
			}
		}
		CHECK(found);
	}
}

TEST_CASE("kuhn_poker_4") {
	auto result = equilibrium::compute<kuhn_poker::State<4>>();
	REQUIRE(result.player_strategies.size() == 2);

	auto result2 = equilibrium_sequential::compute<kuhn_poker::State<4>>();
	CHECK(Approx(result.value) == result2.value);
}

TEST_CASE("stengel") {
	auto result = equilibrium::compute<stengel::State>();
	CHECK(Approx(result.value) == 13);

	auto result_sequentual = equilibrium_sequential::compute<stengel::State>();
	CHECK(Approx(result_sequentual.value) == 13);
}

TEST_CASE("stengel_matrices") {
	// Checks that the same matrices are obtained as in the paper
	// Fast Algorithms for Finding Randomized Strategies in Game Trees
	// http://www.maths.lse.ac.uk/personal/stengel/TEXTE/stoc.pdf
	//
	// (up to permutations)

	equilibrium_sequential::RecurseInfo<stengel::State> recurse_info;
	equilibrium_sequential::all_sequences(&recurse_info);

	Eigen::SparseMatrix<double> E;
	vector<double> e;
	Eigen::SparseMatrix<double> F;
	vector<double> f;
	compute_matrices(recurse_info, &E, &e, &F, &f);
	Eigen::SparseMatrix<double> A;
	compute_matrix_A(recurse_info, &A);

	// Matrices from the paper.
	// clang-format off
	Eigen::MatrixXd E_expected(4, 7);
	E_expected <<
	     1, 0, 0, 0, 0, 0, 0,
		-1, 1, 1, 0, 0, 0, 0,
		0, 0, -1, 1, 1, 0, 0,
		-1, 0, 0, 0, 0, 1, 1;
	Eigen::MatrixXd F_expected(3, 5);
	F_expected << 1, 0, 0, 0, 0,
		-1, 1, 1, 0, 0,
		-1, 0, 0, 1, 1;
	Eigen::MatrixXd A_expected(7, 5);
	A_expected <<
	    2, 0, 0, 0, 0,
		1, 0, 0, 0, 0,
		0, 0, 0, 0, 0,
		0, 2, 4, 0, 0,
		0, 3, -1, 0, 0,
		2, 0, 0, 4, 10,
		4, 0, 0, 6, 3;
	// clang-format on

	CHECK(A_expected.size() == A.size());
	CHECK(A_expected.trace() == A.toDense().trace());
	CHECK(A_expected.sum() == A.sum());
	CHECK(E_expected.size() == E.size());
	CHECK((E_expected.transpose() * E_expected).trace() == (E.toDense().transpose() * E).trace());
	CHECK(E_expected.sum() == E.sum());
	CHECK(F_expected.size() == F.size());
	CHECK((F_expected.transpose() * F_expected).trace() == (F.toDense().transpose() * F).trace());
	CHECK(F_expected.sum() == F.sum());
}

TEST_CASE("leduc_holdem") {
	auto result = equilibrium_sequential::compute<leduc_holdem::State>();
	// Assertion from the paper “Heads-up Limit Hold’em Poker is Solved.”
	CHECK((result.player_strategies[0].size() + result.player_strategies[0].size()) == 288);
}
