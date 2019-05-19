#pragma once

#include <array>
#include <random>
#include <utility>

namespace minimum {
namespace core {
// Initializes an std::mt19937[_64] distribution.
// https://codereview.stackexchange.com/questions/114066/simple-random-generator
// http://www.pcg-random.org/posts/cpp-seeding-surprises.html
template <class Engine>
Engine seeded_engine() {
	std::array<typename Engine::result_type, Engine::state_size> seed_data;
	std::random_device source;
	std::generate(std::begin(seed_data), std::end(seed_data), std::ref(source));
	std::seed_seq seeds(std::begin(seed_data), std::end(seed_data));
	Engine engine(seeds);
	return engine;
}

template <class Engine, typename... Ts>
Engine repeatably_seeded_engine(Ts... is) {
	std::seed_seq seeds = {is...};
	Engine engine(seeds);
	return engine;
}
}  // namespace core
}  // namespace minimum
