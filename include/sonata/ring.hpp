#pragma once
#include "field.hpp"
#include <cstdlib>
#include <functional>
#include <memory>
namespace sonata {
using Ring = std::array<U, 64>;
using RV = std::vector<Ring>;
using SV = std::vector<U>;
using I = int64_t;
using IW = __int128_t;
inline constexpr size_t degree = 64, projrows = 384;
inline constexpr U iv2 = (q + 1) / 2;
inline U qm(I x) { return x >= 0 ? U(x) % q : (q - U(-x) % q) % q; }
inline I center(U x) { return x > q / 2 ? I(x) - I(q) : I(x); }
inline U redwide(W x) {
  W y = (x & mask) + 59 * (x >> 48);
  y = (y & mask) + 59 * (y >> 48);
  U z = U(y);
  z = (z & mask) + 59 * (z >> 48);
  z = (z & mask) + 59 * (z >> 48);
  return z >= q ? z - q : z;
}
inline Ring ra(Ring x, const Ring &y) {
  for (int t = 0; t < 64; t++)
    x[t] = add(x[t], y[t]);
  return x;
}
inline Ring rs(Ring x, const Ring &y) {
  for (int t = 0; t < 64; t++)
    x[t] = sub(x[t], y[t]);
  return x;
}
inline Ring scaleR(Ring x, U a) {
  for (auto &z : x)
    z = mul(z, a);
  return x;
}
inline bool rzero(const Ring &r) {
  return std::all_of(r.begin(), r.end(), [](U x) { return !x; });
}
inline Ring rm(const Ring &a, const Ring &b) {
  W z[127]{};
  for (int i = 0; i < 64; i++)
    if (a[i])
      for (int j = 0; j < 64; j++)
        if (b[j])
          z[i + j] += W(a[i]) * b[j];
  Ring out{};
  for (int i = 0; i < 64; i++)
    out[i] = sub(redwide(z[i]), i + 64 < 127 ? redwide(z[i + 64]) : 0);
  return out;
}
inline Ring dotR(const RV &a, const RV &b) {
  check(a.size() == b.size(), "ring dot dimensions");
  W z[127]{};
  for (size_t k = 0; k < a.size(); k++)
    for (int i = 0; i < 64; i++)
      if (a[k][i])
        for (int j = 0; j < 64; j++)
          if (b[k][j])
            z[i + j] += W(a[k][i]) * b[k][j];
  Ring r{};
  for (int i = 0; i < 64; i++)
    r[i] = sub(redwide(z[i]), i + 64 < 127 ? redwide(z[i + 64]) : 0);
  return r;
}
inline Ring dual(const Ring &x) {
  Ring y{};
  y[0] = x[0];
  for (int t = 1; t < 64; t++)
    y[64 - t] = x[t] ? q - x[t] : 0;
  return y;
}
inline U ctDot(const RV &a, const RV &b) {
  check(a.size() == b.size(), "ct dot length");
  U out = 0;
  for (size_t j = 0; j < a.size(); j++) {
    out = add(out, mul(a[j][0], b[j][0]));
    for (int t = 1; t < 64; t++)
      out = sub(out, mul(a[j][t], b[j][64 - t]));
  }
  return out;
}
inline U cv(const RV &x, size_t i) {
  check(i / 64 < x.size(), "coefficient access");
  return x[i / 64][i % 64];
}
inline void put(RV &x, size_t i, U v) {
  check(i / 64 < x.size(), "coefficient put");
  x[i / 64][i % 64] = v;
}
inline W energy(const RV &x) {
  W e = 0;
  for (const auto &r : x)
    for (U a : r) {
      I z = center(a);
      e += W(z < 0 ? -z : z) * (z < 0 ? -z : z);
    }
  return e;
}
inline RV slice(const RV &v, size_t p, size_t n) {
  check(p + n <= v.size(), "slice bounds");
  return RV(v.begin() + p, v.begin() + p + n);
}
inline void append(SV &out, const RV &v) {
  for (auto &r : v)
    out.insert(out.end(), r.begin(), r.end());
}
inline RV chunks(const SV &v, size_t p, size_t n, bool lift = false) {
  check(p + 64 * n <= v.size(), "chunk bounds");
  RV r(n);
  for (size_t i = 0; i < n; i++) {
    std::copy_n(v.begin() + p + 64 * i, 64, r[i].begin());
    if (lift)
      r[i] = dual(r[i]);
  }
  return r;
}
// Polynomial extended Euclid; inverses are used only modulo q, never in norm
// estimates.
using Poly = SV;
inline void trim(Poly &a) {
  while (!a.empty() && !a.back())
    a.pop_back();
}
inline std::pair<Poly, Poly> divideP(Poly a, Poly b) {
  trim(a);
  trim(b);
  check(!b.empty(), "polynomial zero divisor");
  Poly d(a.size() >= b.size() ? a.size() - b.size() + 1 : 0);
  U invlead = powq(b.back(), q - 2);
  while (a.size() >= b.size()) {
    size_t s = a.size() - b.size();
    U c = mul(a.back(), invlead);
    d[s] = c;
    for (size_t j = 0; j < b.size(); j++)
      a[s + j] = sub(a[s + j], mul(c, b[j]));
    trim(a);
  }
  return {d, a};
}
inline Poly pmul(const Poly &a, const Poly &b) {
  if (a.empty() || b.empty())
    return {};
  Poly c(a.size() + b.size() - 1);
  for (size_t i = 0; i < a.size(); i++)
    for (size_t j = 0; j < b.size(); j++)
      c[i + j] = add(c[i + j], mul(a[i], b[j]));
  trim(c);
  return c;
}
inline Poly psub(Poly a, const Poly &b) {
  a.resize(std::max(a.size(), b.size()));
  for (size_t i = 0; i < b.size(); i++)
    a[i] = sub(a[i], b[i]);
  trim(a);
  return a;
}
inline Ring inverseR(const Ring &a) {
  Poly r0(65);
  r0[0] = r0[64] = 1;
  Poly r1(a.begin(), a.end());
  trim(r1);
  Poly s0{}, s1{1};
  while (!r1.empty()) {
    auto [d, r] = divideP(r0, r1);
    auto s = psub(s0, pmul(d, s1));
    r0 = r1;
    r1 = r;
    s0 = s1;
    s1 = s;
  }
  check(r0.size() == 1 && r0[0], "nonunit challenge");
  U u = powq(r0[0], q - 2);
  Ring out{};
  for (size_t i = 0; i < s0.size(); i++) {
    U z = mul(s0[i], u);
    out[i % 64] = (i / 64) % 2 ? sub(out[i % 64], z) : add(out[i % 64], z);
  }
  Ring oneR{};
  oneR[0] = 1;
  check(rm(a, out) == oneR, "inverse verification");
  return out;
}
inline RV digits(const RV &v, int rho, int ell) {
  RV out(v.size() * ell);
  for (size_t j = 0; j < v.size(); j++)
    for (int t = 0; t < 64; t++) {
      I x = center(v[j][t]), sgn = x < 0 ? -1 : 1;
      x = std::llabs(x);
      for (int l = 0; l < ell; l++) {
        I nx = (x + rho / 2 - 1) / rho;
        out[j * ell + l][t] = qm(sgn * (x - rho * nx));
        x = nx;
      }
      check(x == 0, "digit overflow");
    }
  return out;
}
inline RV recompose(const RV &v, int rho, int ell) {
  check(v.size() % ell == 0, "digit segment");
  RV out(v.size() / ell);
  for (size_t j = 0; j < out.size(); j++) {
    U p = 1;
    for (int l = 0; l < ell; l++) {
      out[j] = ra(out[j], scaleR(v[j * ell + l], p));
      p = mul(p, rho);
    }
  }
  return out;
}
struct Config {
  size_t s, n, rho, raw, N;
  U S, G;
  int b, h;
};
// All roles are column-major prefix views of the same supplied CRS.
struct Master {
  RV data;
  const Ring &at(size_t rows, size_t r, size_t c) const {
    check(rows > 0 && r < rows && c < data.size() / rows, "CRS prefix view");
    return data[c * rows + r];
  }
};
struct Key {
  std::shared_ptr<const Master> master;
  size_t rows, cols;
  RV apply(const RV &x) const {
    check(x.size() == cols, "key input");
    RV y(rows);
    for (size_t r = 0; r < rows; r++) {
      RV a(cols);
      for (size_t c = 0; c < cols; c++)
        a[c] = master->at(rows, r, c);
      y[r] = dotR(a, x);
    }
    return y;
  }
  RV transpose(const RV &w) const {
    check(w.size() == rows, "key transposition");
    RV z(cols);
    for (size_t c = 0; c < cols; c++)
      for (size_t r = 0; r < rows; r++)
        z[c] = ra(z[c], rm(w[r], master->at(rows, r, c)));
    return z;
  }
};
struct System {
  size_t raw;
  SV y;
  std::function<SV(const RV &)> apply;
  std::function<RV(const SV &)> adj;
};

} // namespace sonata
