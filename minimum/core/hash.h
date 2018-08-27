#pragma once

#include <functional>
#include <type_traits>

namespace minimum {
namespace core {

namespace {
std::size_t hash_combine(std::size_t h1, std::size_t h2) {
	// From boost hash_combine.
	h2 ^= h1 + 0x9e3779b9 + (h2 << 6) + (h2 >> 2);
	return h2;
}
}  // namespace

namespace {
// hasher(a, b, c) hashes all its arguments into a single hash.
inline std::size_t hasher() { return 0; }
}  // namespace

// hasher(a, b, c) hashes all its arguments into a single hash.
template <typename T, typename... Ts>
std::size_t hasher(const T& t, const Ts&... ts) {
	std::hash<std::remove_cv_t<std::remove_reference_t<T>>> first_hash;
	return hash_combine(first_hash(t), hasher(ts...));
}
}  // namespace core
}  // namespace minimum

namespace std {
template <typename T, typename U>
struct hash<pair<T, U>> {
	hash<T> hasht;
	hash<U> hashu;
	size_t operator()(const pair<T, U>& p) const {
		return minimum::core::hash_combine(hasht(p.first), hashu(p.second));
	}
};
}  // namespace std
