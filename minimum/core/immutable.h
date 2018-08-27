#pragma once
// See https://github.com/patmorin/arraylayout
// for many more options.

#include <memory>

namespace minimum {
namespace core {

namespace internal {
template <typename T, typename Index>
inline Index find(const T* data, const T& value, const Index n) {
	Index i = 0;
	while (i < n) {
		if (value < data[i]) {
			// Go to left child.
			i = 2 * i + 1;
		} else if (value > data[i]) {
			// Go to right child.
			i = 2 * i + 2;
		} else {
			return i;
		}
	}
	return n;
}
}  // namespace internal

template <typename T, typename Index = int>
class ImmutableSet {
   public:
	using iterator = const T*;
	using const_iterator = const T*;

	// Constructs an ImmutableSet by taking a copy of data.
	//
	// itr: Iterator to the first element of *sorted* data.
	// n: Number of elements in the collection.
	template <typename ForwardIterator>
	ImmutableSet(ForwardIterator itr, Index n_) : n(n_), data(std::make_unique<T[]>(n_)) {
		copy(itr, 0);
	}

	Index size() const { return n; }

	iterator begin() const { return data.get(); }

	iterator end() const { return data.get() + n; }

	iterator find(const T& value) const { return begin() + internal::find(data.get(), value, n); }

	bool contains(const T& value) const { return internal::find(data.get(), value, n) != n; }

	int count(const T& value) const { return contains(value) ? 1 : 0; }

   private:
	template <typename ForwardIterator>
	ForwardIterator copy(ForwardIterator itr, Index i) {
		// https://github.com/patmorin/arraylayout/blob/master/src/eytzinger_array.h
		if (i >= n) {
			return itr;
		}
		// Visit left child
		itr = copy(itr, 2 * i + 1);
		// Put data at the root
		data[i] = *itr;
		itr++;
		// Visit right child
		itr = copy(itr, 2 * i + 2);
		return itr;
	}

	const Index n;
	std::unique_ptr<T[]> data;
};

template <typename K, typename V, typename Index = int>
class ImmutableMap {
   public:
	class iterator {
	   public:
		iterator(const K* keys_, const V* values_) : keys(keys_), values(values_) {}

		std::pair<const K&, const V&> operator*() const { return {*keys, *values}; }

		bool operator!=(const iterator& rhs) const { return keys != rhs.keys; }

		bool operator==(const iterator& rhs) const { return keys == rhs.keys; }

		void operator++() {
			keys++;
			values++;
		}

		void operator+=(Index index) {
			keys += index;
			values += index;
		}

		struct ArrowHelper {
			std::pair<const K&, const V&> value;
			ArrowHelper(const K& k, const V& v) : value(k, v) {}
			const std::pair<const K&, const V&>* operator->() const { return &value; }
		};
		ArrowHelper operator->() const { return ArrowHelper{*keys, *values}; }

	   private:
		const K* keys;
		const V* values;
	};
	using const_iterator = const iterator;

	// Constructs an ImmutableMap by taking a copy of data.
	//
	// itr: Iterator to the first element of *sorted* data.
	//      For example obtained by std::map<K, V>::begin().
	// n: Number of elements in the collection.
	template <typename ForwardIterator>
	ImmutableMap(ForwardIterator itr, Index n_)
	    : n(n_), keys(std::make_unique<K[]>(n_)), values(std::make_unique<V[]>(n_)) {
		copy(itr, 0);
	}

	Index size() const { return n; }

	iterator begin() const { return {keys.get(), values.get()}; }

	iterator end() const { return {keys.get() + n, values.get() + n}; }

	iterator find(const K& key) const {
		iterator itr = begin();
		itr += internal::find(keys.get(), key, n);
		return itr;
	}

	const V& operator[](const K& key) const {
		auto itr = find(key);
		if (itr == end()) {
			throw std::out_of_range("Missing key in ImmutableMap.");
		}
		return (*itr).second;
	}

	const V& at(const K& key) const { return operator[](key); }

	int count(const K& key) const { return internal::find(keys.get(), key, n) != n ? 1 : 0; }

   private:
	template <typename ForwardIterator>
	ForwardIterator copy(ForwardIterator itr, Index i) {
		// https://github.com/patmorin/arraylayout/blob/master/src/eytzinger_array.h
		if (i >= n) {
			return itr;
		}
		// Visit left child
		itr = copy(itr, 2 * i + 1);
		// Put data at the root
		keys[i] = itr->first;
		values[i] = itr->second;
		itr++;
		// Visit right child
		itr = copy(itr, 2 * i + 2);
		return itr;
	}

	const Index n;
	std::unique_ptr<K[]> keys;
	std::unique_ptr<V[]> values;
};
}  // namespace core
}  // namespace minimum