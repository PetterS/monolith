#pragma once

#include <algorithm>
#include <array>
#include <initializer_list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace google {
namespace protobuf {
template <typename T>
class RepeatedField;
}
}  // namespace google

namespace minimum {
namespace core {

//
// Like Python's " ".join([1, 2, 3])
//
template <typename Container>
std::string join(std::string_view joiner, const Container& container);
template <typename Container>
std::string join(char ch, const Container& container);
template <typename T>
std::string join(char ch, const std::initializer_list<T>& container);
template <typename T>
std::string join(std::string_view joiner, const std::initializer_list<T>& container);

namespace {
inline bool contains(std::string_view str, std::string_view substr) {
	return str.find(substr) != std::string_view::npos;
}
}  // namespace

// to_string converts all its arguments to a string and
// concatenates.
//
// Anonymous namespace is needed because the base case for the
// template recursion is a normal function.
//
// Works for
//
//  std::pair
//  std::tuple
//  std::vector
//  std::set
//  std::map
//
// and combinations thereof, e.g. vector<pair<int, string>>.
template <typename... Args>
std::ostream& operator<<(std::ostream& stream, const std::tuple<Args...>& t);

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& stream, const std::pair<T1, T2>& p) {
	stream << '(' << p.first << ", " << p.second << ')';
	return stream;
}

template <typename Tuple, std::size_t i = 0>
void print_tuple(std::ostream& stream, const Tuple& t) {
	if constexpr (i < std::tuple_size<Tuple>::value - 1) {
		stream << std::get<i>(t) << ", ";
		print_tuple<Tuple, i + 1>(stream, t);
	} else {
		stream << std::get<std::tuple_size<Tuple>::value - 1>(t);
	}
}

template <typename... Args>
std::ostream& operator<<(std::ostream& stream, const std::tuple<Args...>& t) {
	stream << '(';
	print_tuple(stream, t);
	stream << ')';
	return stream;
}

template <typename T, typename Alloc>
std::ostream& operator<<(std::ostream& stream, const std::vector<T, Alloc>& v) {
	stream << '[' << join(", ", v) << ']';
	return stream;
}

template <typename T, std::size_t sz>
std::ostream& operator<<(std::ostream& stream, const std::array<T, sz>& v) {
	stream << '[' << join(", ", v) << ']';
	return stream;
}

template <typename T, typename Comp, typename Alloc>
std::ostream& operator<<(std::ostream& stream, const std::set<T, Comp, Alloc>& s) {
	stream << '{' << join(", ", s) << '}';
	return stream;
}

template <typename T, typename Hash, typename Comp, typename Alloc>
std::ostream& operator<<(std::ostream& stream, const std::unordered_set<T, Hash, Comp, Alloc>& s) {
	stream << '{' << join(", ", s) << '}';
	return stream;
}

template <typename T1, typename T2, typename Comp, typename Alloc>
std::ostream& operator<<(std::ostream& stream, const std::map<T1, T2, Comp, Alloc>& m) {
	stream << '[' << join(", ", m) << ']';
	return stream;
}

template <typename T1, typename T2, typename Hash, typename Comp, typename Alloc>
std::ostream& operator<<(std::ostream& stream,
                         const std::unordered_map<T1, T2, Hash, Comp, Alloc>& m) {
	stream << '[' << join(", ", m) << ']';
	return stream;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const google::protobuf::RepeatedField<T>& m) {
	stream << '[' << join(", ", m) << ']';
	return stream;
}

inline void add_to_stream(std::ostream*) {}

template <typename T, typename... Args>
void add_to_stream(std::ostream* stream, T&& t, Args&&... args) {
	(*stream) << std::forward<T>(t);
	add_to_stream(stream, std::forward<Args>(args)...);
}

namespace {
inline std::string to_string() { return {}; }
}  // namespace

// Overload for string literals.
template <size_t n>
std::string to_string(const char (&c_str)[n]) {
	return {c_str};
}

template <typename... Args>
std::string to_string(Args&&... args) {
	std::ostringstream stream;
	add_to_stream(&stream, std::forward<Args>(args)...);
	return stream.str();
}

template <typename T>
T from_string(const std::string& input_string) {
	std::istringstream input_stream(input_string);
	T t;
	input_stream >> t;
	if (!input_stream) {
		std::ostringstream error;
		error << "Could not parse " << typeid(T).name() << " from \"" << input_string << "\".";
		throw std::runtime_error(error.str());
	}
	return t;
}

template <typename T>
T from_string(const std::string& input_string, T default_value) {
	std::istringstream input_stream(input_string);
	T t;
	input_stream >> t;
	if (!input_stream) {
		t = default_value;
	}
	return t;
}

template <typename Container>
std::string join_helper(std::string_view joiner, const Container& container) {
	std::string output = "";
	bool first = true;
	for (const auto& elem : container) {
		if (!first) {
			output += joiner;
		}
		first = false;
		output += to_string(elem);
	}
	return output;
}

template <typename Container>
std::string join(std::string_view joiner, const Container& container) {
	return join_helper(joiner, container);
}

template <typename Container>
std::string join(char ch, const Container& container) {
	std::string joiner;
	joiner += ch;
	return join_helper(joiner, container);
}

template <typename T>
std::string join(char ch, const std::initializer_list<T>& container) {
	std::string joiner;
	joiner += ch;
	return join_helper(joiner, container);
}

template <typename T>
std::string join(std::string_view joiner, const std::initializer_list<T>& container) {
	return join_helper(joiner, container);
}

template <typename Int>
std::string to_string_with_separator(Int input) {
	if (input < 0) {
		return "-" + to_string_with_separator(-input);
	}

	std::string s;
	auto num = to_string(input);
	std::reverse(begin(num), end(num));
	for (std::size_t i = 0; i < num.size(); ++i) {
		s += num[i];
		if (i % 3 == 2 && i < num.size() - 1) {
			s += ',';
		}
	}
	std::reverse(begin(s), end(s));
	return s;
}

namespace {
inline std::vector<std::string>& split(const std::string& s,
                                       char delim,
                                       std::vector<std::string>& parts) {
	std::istringstream sin(s);
	std::string part;
	while (getline(sin, part, delim)) {
		parts.push_back(part);
	}
	return parts;
}

inline std::vector<std::string> split(const std::string& s, char delim) {
	std::vector<std::string> parts;
	split(s, delim, parts);
	return parts;
}

// Removes all spaces in-place.
inline void remove_spaces(std::string* in, int (*pred)(int) = ::isspace) {
	in->erase(std::remove_if(in->begin(), in->end(), pred), in->end());
}

// Removes spaces before and after in-place.
inline std::string_view strip(std::string_view str, int (*pred)(int) = ::isspace) {
	std::string_view::size_type read = 0;
	// Skip initial space.
	for (; read < str.size() && pred(str[read]); ++read) {
	}
	str.remove_prefix(read);

	std::string::size_type ok_size = 0;
	for (read = 0; read < str.size(); ++read) {
		if (!pred(str[read])) {
			ok_size = read + 1;
		}
	}
	// Resize to skip the last spaces.
	str.remove_suffix(str.size() - ok_size);

	return str;
}

inline bool ends_with(std::string_view s, std::string_view ending) {
	if (s.length() >= ending.length()) {
		return std::equal(s.end() - ending.length(), s.end(), ending.begin());
	} else {
		return false;
	}
}
}  // namespace
}  // namespace core
}  // namespace minimum
