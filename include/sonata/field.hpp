#pragma once
#include "normal_form.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>
namespace sonata {
using U = uint64_t;
using W = __uint128_t;
using E = std::array<U, 4>;
inline constexpr U q = (1ULL << 48) - 59, mask = (1ULL << 48) - 1;
inline U add(U a, U b) {
  U c = a + b;
  return c >= q ? c - q : c;
}
inline U sub(U a, U b) { return a >= b ? a - b : a + q - b; }
inline U mul(U a, U b) {
  W x = W(a) * b;
  U r = U(x & mask) + 59 * U(x >> 48);
  r = (r & mask) + 59 * (r >> 48);
  return r >= q ? r - q : r;
}
inline U powq(U a, U n) {
  U b = 1;
  for (; n; n >>= 1, a = mul(a, a))
    if (n & 1)
      b = mul(b, a);
  return b;
}
inline constexpr E zero{}, one{1, 0, 0, 0};
inline const std::array<U, 42> high_inv = []() {
  std::array<U, 42> a{};
  for (int i = 0; i < 42; i++)
    a[i] = powq(9 + i, q - 2);
  return a;
}();
inline E ea(E a, E b) {
  for (int i = 0; i < 4; i++)
    a[i] = add(a[i], b[i]);
  return a;
}
inline E es(E a, E b) {
  for (int i = 0; i < 4; i++)
    a[i] = sub(a[i], b[i]);
  return a;
}
inline E sc(E a, U b) {
  for (auto &v : a)
    v = mul(v, b);
  return a;
}
inline E red(U *v) {
  return {sub(sub(v[0], mul(2, v[4])), mul(8, v[6])), sub(v[1], mul(2, v[5])),
          add(add(v[2], mul(4, v[4])), mul(14, v[6])), add(v[3], mul(4, v[5]))};
}
inline E em(E a, E b) {
  U v[7]{};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      v[i + j] = add(v[i + j], mul(a[i], b[j]));
  return red(v);
}
inline E sq(E a) {
  U v[7]{};
  for (int i = 0; i < 4; i++) {
    v[2 * i] = add(v[2 * i], mul(a[i], a[i]));
    for (int j = i + 1; j < 4; j++)
      v[i + j] = add(v[i + j], mul(2, mul(a[i], a[j])));
  }
  return red(v);
}
inline void check(bool b, const char *msg) {
  if (!b)
    throw std::runtime_error(msg);
}
struct Lin {
  E w;
  int s, c;
};
struct Transcript {
  std::vector<E> g0, g2, ch;
  E terminal{};
  U mults = 0, squares = 0, foldmults = 0;
};
inline bool verify(const Transcript &t, E claim, E bound_terminal) {
  if (t.g0.size() != t.g2.size() || t.g0.size() != t.ch.size())
    return false;
  for (size_t j = 0; j < t.g0.size(); j++) {
    E g1 = es(es(claim, sc(t.g0[j], 2)), t.g2[j]);
    claim = ea(t.g0[j], em(t.ch[j], ea(g1, em(t.ch[j], t.g2[j]))));
  }
  return bound_terminal == t.terminal && sq(bound_terminal) == claim;
}
inline std::vector<E> eqweights(const std::vector<E> &r) {
  std::vector<E> w(1, one);
  for (auto x : r) {
    size_t n = w.size();
    w.resize(2 * n);
    for (size_t i = 0; i < n; i++) {
      E old = w[i];
      w[i] = em(old, es(one, x));
      w[i + n] = em(old, x);
    }
  }
  return w;
}

} // namespace sonata
