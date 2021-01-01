// Investigates the two-player version of Svälta Räv (War).
//
// Specifically, how often the game can end in a draw, which happened the second
// or third game after I taught the game to my daughter.
//
// Quite fast Monte-Carlo implementation, with zero allocations per played game.
// 2µs per complete game.
//
// Petter Strandmark
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include <absl/random/random.h>

#include <minimum/core/check.h>
#include <minimum/core/color.h>
using namespace std;

enum class Suit { SPADE, CLUB, HEART, DIAMOND };
typedef pair<Suit, int> Card;

// Not needed in this program really, but perhaps useful in the future?
ostream& operator<<(ostream& out, const Card& card) {
	if (card.second == 14) {
		out << "A";
	} else if (card.second == 13) {
		cout << "K";
	} else if (card.second == 12) {
		cout << "Q";
	} else if (card.second == 11) {
		cout << "J";
	} else {
		out << card.second;
	}

	if (card.first == Suit::SPADE) {
		out << "♠";
	} else if (card.first == Suit::CLUB) {
		out << "♣";
	} else if (card.first == Suit::HEART) {
		out << minimum::core::RED << "♥" << minimum::core::NORMAL;
	} else if (card.first == Suit::DIAMOND) {
		out << minimum::core::RED << "♦" << minimum::core::NORMAL;
	}
	return out;
}

int main(int argc, char* argv[]) {
	absl::BitGen rng;

	int iterations = 1000;
	if (argc > 1) {
		iterations = atoi(argv[1]);
	}
	cerr << "Running " << iterations << " iterations.\n";

	vector<Card> cards;
	for (int val = 2; val <= 14; ++val) {
		cards.emplace_back(Suit::SPADE, val);
		cards.emplace_back(Suit::CLUB, val);
		cards.emplace_back(Suit::HEART, val);
		cards.emplace_back(Suit::DIAMOND, val);
	}
	for (auto& card : cards) {
		cerr << card << " ";
	}
	cerr << endl;

	vector<Card> player1_draw, player1_spoils, player2_draw, player2_spoils, current;
	player1_draw.reserve(52);
	player1_spoils.reserve(52);
	player2_draw.reserve(52);
	player2_spoils.reserve(52);
	current.reserve(52);

	int player1_wins = 0;
	int player2_wins = 0;
	int draws = 0;
	int longest_game = -1;
	int shortest_game = 1'000'000'000;
	map<int, int> game_lengths;

	minimum::core::Timer timer("Simulating");
	for (int iter = 1; iter <= iterations; ++iter) {
		player1_draw.clear();
		player1_spoils.clear();
		player2_draw.clear();
		player2_spoils.clear();
		current.clear();

		shuffle(begin(cards), end(cards), rng);
		minimum_core_assert(cards.size() == 52);
		for (int i = 0; i < 26; ++i) {
			player1_draw.push_back(cards[2 * i]);
			player2_draw.push_back(cards[2 * i + 1]);
		}

		int turns = 0;
		while (!player1_draw.empty() && !player2_draw.empty()) {
			turns++;

			auto p1 = player1_draw.back();
			player1_draw.pop_back();
			auto p2 = player2_draw.back();
			player2_draw.pop_back();

			current.push_back(p1);
			current.push_back(p2);
			if (p1.first == p2.first) {
				if (p1.second > p2.second) {
					// cerr << "Player 1 got a pile of " << current.size() << endl;
					player1_spoils.insert(end(player1_spoils), begin(current), end(current));
				} else {
					minimum_core_assert(p1.second != p2.second);
					// cerr << "Player 2 got a pile of " << current.size() << endl;
					player2_spoils.insert(end(player2_spoils), begin(current), end(current));
				}
				current.clear();
			}

			if (player1_draw.empty()) {
				shuffle(begin(player1_spoils), end(player1_spoils), rng);
				player1_draw.swap(player1_spoils);
			}
			if (player2_draw.empty()) {
				shuffle(begin(player2_spoils), end(player2_spoils), rng);
				player2_draw.swap(player2_spoils);
			}
		}

		longest_game = max(longest_game, turns);
		shortest_game = min(shortest_game, turns);
		game_lengths[turns]++;

		if (player1_draw.empty() && player2_draw.empty()) {
			draws++;
		} else if (player1_draw.empty()) {
			// cerr << "Player 2 wins with " << current.size() << " remaining." << endl;
			player2_wins++;
		} else {
			// cerr << "Player 1 wins with " << current.size() << " remaining." << endl;
			player1_wins++;
		}
	}
	timer.OK();

	cout << "Player 1 wins: " << player1_wins << endl;
	cout << "Player 2 wins: " << player2_wins << endl;
	cout << "Draws: " << draws << " (" << 100.0 * double(draws) / iterations << "%)" << endl;
	cout << endl;
	cout << "Longest game: " << longest_game << " turns." << endl;
	cout << "Shortest game: " << shortest_game << " turns." << endl;
	int max_bin = -1;
	for (auto& item : game_lengths) {
		max_bin = max(max_bin, item.second);
	}
	constexpr int histogram_width = 50;
	int n_printed = 0;
	for (auto& item : game_lengths) {
		if (++n_printed > 20) {
			cout << "   ..." << endl;
			break;
		}

		cout << setw(6) << right << item.first << " turns: ";
		int size = (histogram_width * item.second) / max_bin;
		for (int i = 0; i < size; ++i) {
			cout << "█";
		}
		for (int i = 0; i < histogram_width - size; ++i) {
			cout << " ";
		}
		cout << ": " << item.second << " (" << 100.0 * double(item.second) / iterations << "%)"
		     << endl;
	}
}
