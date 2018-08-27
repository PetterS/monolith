#pragma once
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace minimum {
namespace core {

namespace internal {
template <typename Int>
struct RangeIterator {
	Int value;
	bool operator!=(const RangeIterator& rhs) const { return value != rhs.value; }

	Int operator*() const { return value; }

	RangeIterator& operator++() {
		value++;
		return *this;
	}
	RangeIterator& operator++(int) {
		value++;
		return *this;
	}
};

template <typename Int>
struct Range {
	using iterator = RangeIterator<Int>;
	Int begin_;
	Int end_;
	iterator begin() const { return {begin_}; }
	iterator end() const { return {end_}; }
};

template <typename Iterator1, typename Iterator2>
struct ZipIterator2 {
	Iterator1 itr1;
	Iterator2 itr2;

	ZipIterator2(Iterator1 itr1_, Iterator2 itr2_)
	    : itr1(std::move(itr1_)), itr2(std::move(itr2_)) {}

	bool operator!=(const ZipIterator2& rhs) const { return itr1 != rhs.itr1 && itr2 != rhs.itr2; }

	auto operator*() const { return std::forward_as_tuple(*itr1, *itr2); }

	ZipIterator2& operator++() {
		itr1++;
		itr2++;
		return *this;
	}
};

template <typename Container1, typename Container2>
struct Zip {
	Container1 container1;
	Container2 container2;
	Zip(Container1 container1_, Container2 container2_)
	    : container1(container1_), container2(container2_) {}

	ZipIterator2<decltype(container1.begin()), decltype(container2.begin())> begin() {
		return {container1.begin(), container2.begin()};
	}

	ZipIterator2<decltype(container1.end()), decltype(container2.end())> end() {
		return {container1.end(), container2.end()};
	}
};
}  // namespace internal

// Allows loops similar to Python syntax:
//
//    for (int i : range(10)) { ... }
//
template <typename Int>
internal::Range<Int> range(Int n) {
	return {0, n};
}
// Allows loops similar to Python syntax:
//
//    for (int i : range(5, 10)) { ... }
//
template <typename Int>
internal::Range<Int> range(Int begin, Int end) {
	return {begin, end};
}

// Allows zip similar to Python syntax:
//
//    for (auto tuple: zip(container1, container2)) { ... }
//
template <typename... T>
internal::Zip<T...> zip(T&&... containers) {
	return {containers...};
}
}  // namespace core
}  // namespace minimum
