#pragma once

#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdarg>
#include <cstring>
#include <span>
#include <string>
#include <string_view>

inline std::string strf(const char* format, ...) {
  char buffer[1024];
  va_list args;
  va_start(args, format);
  vsprintf(buffer, format, args);
  return std::string(buffer);
};

namespace tokens {

using StringView = std::string_view;
using StringViewSpan = std::span<const StringView>;

inline bool startswith(StringView str, StringView prefix) {
  return (str.size() >= prefix.size()) && (std::memcmp(str.data(), prefix.data(), prefix.size()) == 0);
}

inline bool endswith(StringView str, StringView suffix) {
  return (str.size() >= suffix.size()) && (std::memcmp(str.data() + (str.size() - suffix.size()), suffix.data(), suffix.size()) == 0);
}

struct IsHashComment {
  constexpr bool operator()(char c) const noexcept {
    return c == '#';
  }
};

template <typename Predicate>
inline bool skipline_tmpl(StringView str, Predicate&& predicate) noexcept {
  return str.empty() || predicate(str.front());
}

inline bool skipline(StringView str) noexcept { return skipline_tmpl(str, IsHashComment{}); }

/* ------------------------------------------------------------------------- */

// Struct that holds a fixed size array of token.
// Warning: Tokens are std::string_view so the TokenSet does not own the data
template <size_t Size> struct TokenSetTmpl {
  alignas(64) std::array<StringView, Size> tokens{};
  size_t len = 0;

  inline void clear() noexcept { len = 0; }

  inline StringView* data() noexcept { return tokens.data(); }
  inline const StringView* data() const noexcept { return tokens.data(); }

  inline StringView& operator[](size_t i) noexcept { return tokens[i]; }
  inline const StringView& operator[](size_t i) const noexcept { return tokens[i]; }

  template <size_t index> inline constexpr StringView& at() noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  template <size_t index> inline constexpr const StringView& at() const noexcept {
    static_assert(index < Size);
    return tokens[index];
  }

  constexpr size_t size() const noexcept { return len; }
  constexpr size_t capacity() const noexcept { return Size; }
  constexpr inline bool full() const noexcept { return size() >= capacity(); }

  inline void push_back(StringView token) noexcept {
    if (!full()) tokens[len++] = token;
  }

  // iteration accessor
  StringView* begin() noexcept { return tokens.begin(); }
  const StringView* begin() const noexcept { return tokens.begin(); }
  StringView* end() noexcept { return tokens.begin() + len; }
  const StringView* end() const noexcept { return tokens.begin() + len; }
};

using TokenSet = TokenSetTmpl<256>;

/* ------------------------------------------------------------------------- */

struct IsSpace {
  constexpr bool operator()(char c) const noexcept {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
  }
};

struct SpaceDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ' ' || c == '\t'; }
};

struct CommaDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ','; }
};

struct ColumnDelimiter {
  constexpr bool operator()(char c) const noexcept { return c == ':'; }
};

template <typename Predicate>
inline StringView ltrim(StringView str, Predicate&& predicate) noexcept {
  while (!str.empty() && predicate(str.front())) {
    str.remove_prefix(1);
  }
  return str;
}

template <typename Predicate>
inline StringView rtrim(StringView str, Predicate&& predicate) noexcept {
  while (!str.empty() && predicate(str.back())) {
    str.remove_suffix(1);
  }
  return str;
}

template <typename Predicate>
inline StringView trim(StringView str, Predicate&& predicate) noexcept {
  while (!str.empty() && predicate(str.front())) {
    str.remove_prefix(1);
  }
  while (!str.empty() && predicate(str.back())) {
    str.remove_suffix(1);
  }
  return str;
}

inline StringView trim_spaces(StringView str) noexcept { return trim(str, IsSpace{}); }

template <typename Callback, typename Delimiter>
inline void split(StringView str, Delimiter&& delimiter, Callback&& callback) noexcept {
  const char* data = str.data();
  const char* end = data + str.size();
  while (data < end) {
    while (data < end && delimiter(*data)) {
      ++data;
    }
    const char* token_start = data;
    while (data < end && !delimiter(*data)) {
      ++data;
    }
    if (token_start != data) {
      callback(StringView(token_start, data - token_start));
    }
  }
}

template <typename Delimiter, size_t N>
inline void tokenize_tmpl(const StringView str, TokenSetTmpl<N>& tokens, Delimiter&& delimiter) {
  tokens.clear();
  split(str, delimiter, [&tokens](StringView token) {
    if (!token.empty()) {
      tokens.push_back(token);
    }
  });
}

template <typename Delimiter, size_t N>
inline void tokenize_nmax_tmpl(const StringView str, TokenSetTmpl<N>& tokens, Delimiter&& delimiter,
                               size_t nmax) {
  tokens.clear();
  split(str, delimiter, [&tokens, &nmax](StringView token) {
    if (!token.empty() && tokens.size() < nmax) {
      tokens.push_back(token);
    }
  });
}

inline void tokenize(const StringView str, TokenSet& tokens) {
  return tokenize_tmpl(str, tokens, SpaceDelimiter{});
}

inline void tokenize(const StringView str, TokenSet& tokens, size_t nmax) {
  return tokenize_nmax_tmpl(str, tokens, SpaceDelimiter{}, nmax);
}

/* ------------------------------------------------------------------------- */

// Return true if 'token' match any of the needles
// Stop at the first match and set 'index' to the matching index in 'needles'.
inline bool token_match_any(StringView token, StringViewSpan needles, size_t& index) {
  const size_t nc = needles.size();
  if (nc == 0) {
    return false;
  }
  for (size_t i = 0; i < nc; ++i) {
    if (token == needles[i]) {
      index = i;
      return true;
    }
  }
  return false;
}

inline bool token_match_any(StringView token, StringViewSpan needles) {
  size_t index;
  return token_match_any(token, needles, index);
}

// Return true when any token in 'tokens' match any needle in 'needles'
// Stop at the first match and provide i and j indices of the token and needles respectively.
inline bool token_set_match_any(StringViewSpan tokens, StringViewSpan needles, size_t& ii, size_t& jj) {
  const size_t token_count = tokens.size();
  const size_t needle_count = needles.size();

  if (needle_count == 0) {
    return false;
  }

  for (size_t i = 0; i < token_count; ++i) {
    for (size_t j = 0; j < needle_count; ++j) {
      if (tokens[i] == tokens[j]) {
        ii = i; // token index
        jj = j; // needle index
        return true;
      }
    }
  }

  return false;
}

inline bool token_set_match_any(StringViewSpan tokens, StringViewSpan needles) {
  size_t i, j;
  return token_set_match_any(tokens, needles, i, j);
}

// Returns true if the sequence of needles is present in the sequence of tokens
inline bool token_set_match_seq(StringViewSpan tokens, StringViewSpan needles, size_t& index) {
  const size_t token_count = tokens.size();
  const size_t needle_count = needles.size();

  if (needle_count == 0 || token_count < needle_count) {
    return false;
  }

  for (size_t i = 0; i <= token_count - needle_count; ++i) {
    bool match = true;
    for (size_t j = 0; j < needle_count; ++j) {
      if (tokens[i + j] != needles[j]) {
        match = false;
        break;
      }
    }
    if (match) {
      index = i;
      return true;
    }
  }
  return false;
}

inline bool token_set_match_seq(StringViewSpan tokens, StringViewSpan needles) {
  size_t index;
  return token_set_match_seq(tokens, needles, index);
}

template <size_t N> using TokenNeedles = std::array<StringView, N>;

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles) {
  return token_match_any(token, StringViewSpan(needles.data(), N));
}

template <size_t N> inline bool token_match_any(const StringView token, const TokenNeedles<N>& needles, size_t& index) {
  return token_match_any(token, StringViewSpan(needles.data(), N), index);
}

template <size_t N, size_t M> inline bool token_set_match_any(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles) {
  return token_set_match_any(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N));
}

template <size_t N, size_t M> inline bool token_set_match_any(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles, size_t& i, size_t& j) {
  return token_set_match_any(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N), i, j);
}

template <size_t N, size_t M> inline bool token_set_match_seq(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles) {
  return token_set_match_seq(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N));
}

template <size_t N, size_t M> inline bool token_set_match_seq(const TokenSetTmpl<M>& tokens, const TokenNeedles<N>& needles, size_t& index) {
  return token_set_match_seq(StringViewSpan(tokens.data(), tokens.size()), StringViewSpan(needles.data(), N), index);
}

/* ------------------------------------------------------------------------- */
// convertion functions

// convert a string_view token to numeric type
template <typename T> bool token_to_num(const StringView token, T& rval) {
  auto [ptr, err] = std::from_chars(token.begin(), token.begin() + token.size(), rval);
  return (err == std::errc() && ptr == token.end());
}

// Dark magic that expand to tokens_to_num with automatic type infering
#define parb ()
#define IOexpand(...) IOexpandB(IOexpandB(IOexpandB(IOexpandB(__VA_ARGS__))))
#define IOexpandB(...) IOexpandA(IOexpandA(IOexpandA(IOexpandA(__VA_ARGS__))))
#define IOexpandA(...) __VA_ARGS__

#define call_parse(token, value) token_to_num<decltype(value)>(token, value)
#define IOfor_each(macro, ...) __VA_OPT__(IOexpand(IOfor_each_helper(macro, __VA_ARGS__)))
#define IOfor_each_helper(macro, token, value, ...)                                                \
  macro(token, value) __VA_OPT__(IOfor_each_again parb(macro, __VA_ARGS__))
#define IOfor_each_again() &&IOfor_each_helper
#define tokens_to_num(...) IOfor_each(call_parse, __VA_ARGS__)

} // namespace tokens

namespace cmp {

// implementation from c++ reference to be use without c++20

template <class T, class U> constexpr bool cmp_equal_int(T t, U u) noexcept {
  if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
    return t == u;
  } else if constexpr (std::is_signed_v<T>) {
    // T is signed, U is unsigned
    return t >= 0 && static_cast<std::make_unsigned_t<T>>(t) == u;
  } else {
    // T is signed, U is unsigned
    return u >= 0 && t == static_cast<std::make_unsigned_t<U>>(u);
  }
}

template <class T, class U> constexpr bool cmp_not_equal_int(T t, U u) noexcept {
  return !cmp_equal_int(t, u);
}

template <class T, class U> constexpr bool cmp_less_int(T t, U u) noexcept {
  if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
    return t < u;
  } else if constexpr (std::is_signed_v<T>) {
    return t < 0 || static_cast<std::make_unsigned_t<T>>(t) < u;
  } else {
    return u >= 0 && t < static_cast<std::make_unsigned_t<U>>(u);
  }
}

template <class T, class U> constexpr bool cmp_greater_int(T t, U u) noexcept {
  return cmp_less_int(u, t);
}

template <class T, class U> constexpr bool cmp_less_equal_int(T t, U u) noexcept {
  return !cmp_less_int(u, t);
}

template <class T, class U> constexpr bool cmp_greater_equal_int(T t, U u) noexcept {
  return !cmp_less_int(t, u);
}

// wrapper functions

template <class T, class U> constexpr inline bool eq(T t, U u) noexcept {
  return cmp_equal_int(t, u);
}
template <class T, class U> constexpr inline bool ne(T t, U u) noexcept {
  return cmp_not_equal_int(t, u);
}
template <class T, class U> constexpr inline bool lt(T t, U u) noexcept {
  return cmp_less_int(t, u);
}
template <class T, class U> constexpr inline bool le(T t, U u) noexcept {
  return cmp_less_equal_int(t, u);
}
template <class T, class U> constexpr inline bool gt(T t, U u) noexcept {
  return cmp_greater_int(t, u);
}
template <class T, class U> constexpr inline bool ge(T t, U u) noexcept {
  return cmp_greater_equal_int(t, u);
}

template <class T, class U> constexpr auto min(const T& a, const U& b) { return lt(b, a) ? b : a; }

template <class T, class U> constexpr auto max(const T& a, const U& b) { return lt(b, a) ? a : b; }

template <class T, class U> constexpr inline bool nearly_eq(T t, U u,
             typename std::common_type<T, U>::type epsilon =
                 std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {

  using common_t = typename std::common_type<T, U>::type;
  static_assert(std::is_floating_point_v<common_t>,
                "nearly_equal requires at least one floating-point type");

  common_t a = static_cast<common_t>(t);
  common_t b = static_cast<common_t>(u);
  common_t diff = std::fabs(a - b);
  common_t norm = std::max(std::fabs(a), std::fabs(b));
  return diff <= epsilon * norm;
}

template <class T, class U> constexpr inline bool nearly_lt(T t, U u,
             typename std::common_type<T, U>::type epsilon =
                 std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {

  using common_t = typename std::common_type<T, U>::type;
  return !nearly_eq(t, u, epsilon) && static_cast<common_t>(t) < static_cast<common_t>(u);
}

template <class T, class U> constexpr inline bool nearly_gt(T t, U u,
             typename std::common_type<T, U>::type epsilon =
                 std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return nearly_lt(u, t, epsilon);
}

template <class T, class U> constexpr inline bool nearly_le(T t, U u,
             typename std::common_type<T, U>::type epsilon =
                 std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return !nearly_lt(t, u, epsilon);
}

template <class T, class U> constexpr inline bool nearly_ge(T t, U u,
             typename std::common_type<T, U>::type epsilon =
                 std::numeric_limits<typename std::common_type<T, U>::type>::epsilon()) noexcept {
  return !nearly_lt(u, t, epsilon);
}

template <typename T> constexpr T floor(T x) {
  static_assert(std::is_floating_point_v<T>, "constexpr_floor: T must be floating-point");
  long long xi = static_cast<long long>(x); // truncates toward zero
  return (x < static_cast<T>(xi)) ? static_cast<T>(xi - 1) : static_cast<T>(xi);
}

template <typename T> constexpr T wrap(T x, T max) {
  static_assert(std::is_floating_point_v<T>, "constexpr_wrap: T must be floating-point");
  T r = x - max * cmp::floor(x / max);
  return (r < static_cast<T>(0)) ? r + max : r;
}

}
