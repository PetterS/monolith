#pragma once
#include <type_traits>
#include <vector>

namespace minimum {
namespace core {

// make_grid<T>(size1, f) creates a vector<T> where each
// element is initialized with f().
template <typename T, typename Lambda>
auto make_grid(typename std::enable_if<!std::is_integral<Lambda>::value, std::size_t>::type s,
               const Lambda& f) {
	std::vector<T> v;
	v.reserve(s);
	for (std::size_t i = 0; i < s; ++i) {
		v.emplace_back(f());
	}
	return v;
}

// make_grid<T>(size1, size2,  f) creates a vector<vector<T>> where each
// element is initialized with f().
template <typename T, typename Lambda>
auto make_grid(typename std::enable_if<!std::is_integral<Lambda>::value, std::size_t>::type s,
               std::size_t s2,
               const Lambda& f) {
	using U = decltype(make_grid<T>(s2, f));
	std::vector<U> v;
	v.reserve(s);
	for (std::size_t i = 0; i < s; ++i) {
		v.emplace_back(make_grid<T>(s2, f));
	}
	return v;
}

// make_grid<T>(size1, size2, size3, f) creates a vector<vector<vector<T>>> where each
// element is initialized with f().
template <typename T, typename Lambda>
auto make_grid(typename std::enable_if<!std::is_integral<Lambda>::value, std::size_t>::type s,
               std::size_t s2,
               std::size_t s3,
               const Lambda& f) {
	using U = decltype(make_grid<T>(s2, s3, f));
	std::vector<U> v;
	v.reserve(s);
	for (std::size_t i = 0; i < s; ++i) {
		v.emplace_back(make_grid<T>(s2, s3, f));
	}
	return v;
}

// make_grid<T>(size1) creates a vector<T> where each
// element is deafult-initialized.
template <typename T>
auto make_grid(std::size_t s) {
	return std::vector<T>(s);
}

// make_grid<T>(size1, size2) creates a vector<vector<T>> where each
// element is deafult-initialized.
template <typename T>
auto make_grid(std::size_t s, std::size_t s2) {
	using U = decltype(make_grid<T>(s2));
	return std::vector<U>(s, make_grid<T>(s2));
}

// make_grid<T>(size1, size2, size3) creates a vector<vector<vector<T>>> where each
// element is deafult-initialized.
template <typename T>
auto make_grid(std::size_t s, std::size_t s2, std::size_t s3) {
	using U = decltype(make_grid<T>(s2, s3));
	return std::vector<U>(s, make_grid<T>(s2, s3));
}

template <typename T, typename Index>
class Grid3D {
   public:
	Grid3D(Index m_, Index n_, Index o_) : o(o_), on(o_ * n_), data(new T[m_ * n_ * o_]) {}
	Grid3D(const Grid3D& other) = delete;
	~Grid3D() { delete[] data; }
	inline T& operator()(Index i, Index j, Index k) { return data[i * on + j * o + k]; }
	inline const T& operator()(Index i, Index j, Index k) const { return data[i * on + j * o + k]; }

   private:
	Index o, on;
	T* data;
};
}  // namespace core
}  // namespace minimum
