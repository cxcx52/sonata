#pragma once
#include "ring.hpp"
namespace sonata {
using Bytes = std::vector<uint8_t>;
using Bits = std::vector<U>;
inline bool bit(const Bits &v, size_t p) { return (v[p / 64] >> (p % 64)) & 1; }
inline U window(const Bits &a, size_t p) {
  size_t w = p / 64, s = p % 64;
  U z = w < a.size() ? a[w] >> s : 0;
  if (s && w + 1 < a.size())
    z |= a[w + 1] << (64 - s);
  return z;
}
struct Projection {
  size_t N;
  int b, h;
  std::vector<int8_t> pi;
  Bytes seed;
  size_t hashcalls = 0;
};
inline Projection expand_projection(const Config &c, const Bytes &seed) {
  check(c.N > 0 && c.N <= SIZE_MAX / (2 * projrows) && c.b > 0 && c.h >= 0,
        "projection configuration");
  size_t seedbits = size_t(c.b) + size_t(c.h) * (3 * size_t(c.b) - 1);
  check(seed.size() == (seedbits + 7) / 8, "complete seed length");
  if (seedbits % 8)
    check(!(seed.back() >> (seedbits % 8)), "seed high padding");
  size_t at = 0;
  auto take = [&](size_t n) {
    Bits z((n + 63) / 64);
    for (size_t j = 0; j < n; j++, at++)
      z[j / 64] |= U((seed[at / 8] >> (at % 8)) & 1) << (j % 64);
    return z;
  };
  Bits initial = take(c.b);
  std::vector<Bits> a, t;
  for (int j = 0; j < c.h; j++) {
    a.push_back(take(2 * c.b - 1));
    t.push_back(take(c.b));
  }
  check(at == seedbits, "seed layout");
  Projection out{c.N, c.b, c.h, {}, seed, 0};
  out.pi.reserve(projrows * c.N);
  int pending = -1;
  auto hash = [&](int j, const Bits &x) {
    out.hashcalls++;
    Bits y = t[j];
    for (int r = 0; r < c.b; r++) {
      U p = 0;
      for (size_t w = 0; w < x.size(); w++)
        p ^= window(a[j], r + 64 * w) & x[w];
      if (__builtin_parityll(p))
        y[r / 64] ^= U(1) << (r % 64);
    }
    return y;
  };
  std::function<void(int, const Bits &)> emit = [&](int h, const Bits &x) {
    if (out.pi.size() == projrows * c.N)
      return;
    if (!h) {
      for (int j = 0; j < c.b && out.pi.size() < projrows * c.N; j++) {
        int v = bit(x, j);
        if (pending < 0)
          pending = v;
        else {
          out.pi.push_back(pending - v);
          pending = -1;
        }
      }
      return;
    }
    emit(h - 1, x);
    if (out.pi.size() < projrows * c.N) {
      auto y = hash(h - 1, x);
      emit(h - 1, y);
    }
  };
  emit(c.h, initial);
  check(out.pi.size() == projrows * c.N && pending < 0,
        "projection expansion length");
  return out;
}
inline SV project(const Projection &p, const RV &w) {
  check(w.size() * 64 == p.N && p.pi.size() == p.N * projrows,
        "projection dimensions");
  SV y(projrows);
  for (size_t i = 0; i < projrows; i++) {
    IW a = 0;
    for (size_t j = 0; j < p.N; j++)
      a += IW(p.pi[i * p.N + j]) * center(cv(w, j));
    I v = I(a % IW(q));
    y[i] = qm(v);
  }
  return y;
}
inline bool proj_norm(const SV &p, U S) {
  if (p.size() != projrows)
    return false;
  W e = 0;
  for (U a : p) {
    I z = center(a);
    e += W(z < 0 ? -z : z) * (z < 0 ? -z : z);
  }
  return e <= W(projrows) * S;
}
inline System align(const System &src, const Config &c) {
  check(src.raw == c.raw && c.raw <= c.N && c.N % 64 == 0,
        "input raw dimensions");
  System dst;
  dst.raw = c.N;
  dst.y = src.y;
  dst.y.resize(src.y.size() + c.N - src.raw);
  dst.apply = [src, c](const RV &w) {
    check(w.size() == c.N / 64, "aligned witness size");
    auto y = src.apply(w);
    for (size_t j = src.raw; j < c.N; j++)
      y.push_back(cv(w, j));
    return y;
  };
  dst.adj = [src, c](const SV &a) {
    check(a.size() == src.y.size() + c.N - src.raw, "aligned row count");
    SV b(a.begin(), a.begin() + src.y.size());
    RV w = src.adj(b);
    w.resize(c.N / 64);
    for (size_t j = src.raw; j < c.N; j++) {
      Ring t{};
      t[j % 64] = a[src.y.size() + j - src.raw];
      w[j / 64] = ra(w[j / 64], dual(t));
    }
    return w;
  };
  return dst;
}
struct Aggregated {
  std::array<RV, 4> phi;
  std::array<U, 4> rhs{};
};
inline Aggregated aggregate(const System &sys, const Projection &pr,
                            const SV &p, E zeta, E chi) {
  check(chi == one, "full-field aggregation has chi=1");
  check(p.size() == projrows && pr.N == sys.raw &&
            pr.pi.size() == projrows * pr.N,
        "projection shape");
  size_t rows = sys.y.size() + projrows;
  std::array<SV, 4> weights;
  for (auto &w : weights)
    w.resize(rows);
  E power = one;
  for (size_t j = 0; j < rows; j += 4) {
    for (size_t t = 0; t < 4 && j + t < rows; t++) {
      E basis{};
      basis[t] = 1;
      auto v = em(chi, em(power, basis));
      for (int e = 0; e < 4; e++)
        weights[e][j + t] = v[e];
    }
    power = em(power, zeta);
  }
  Aggregated ag;
  for (int e = 0; e < 4; e++) {
    auto &a = weights[e];
    SV aa(a.begin(), a.begin() + sys.y.size());
    ag.phi[e] = sys.adj(aa);
    ag.phi[e].resize(pr.N / 64);
    SV coeff(pr.N);
    for (size_t r = 0; r < projrows; r++) {
      U weight = a[sys.y.size() + r];
      ag.rhs[e] = add(ag.rhs[e], mul(weight, p[r]));
      const int8_t *row = pr.pi.data() + r * pr.N;
      for (size_t j = 0; j < pr.N; j++) {
        if (row[j] > 0)
          coeff[j] = add(coeff[j], weight);
        else if (row[j] < 0)
          coeff[j] = sub(coeff[j], weight);
      }
    }
    for (size_t j = 0; j < pr.N / 64; j++) {
      Ring v{};
      std::copy_n(coeff.begin() + 64 * j, 64, v.begin());
      ag.phi[e][j] = ra(ag.phi[e][j], dual(v));
    }
    for (size_t j = 0; j < sys.y.size(); j++)
      ag.rhs[e] = add(ag.rhs[e], mul(a[j], sys.y[j]));
  }
  return ag;
}
inline size_t tri(size_t s, size_t a, size_t b) {
  check(a <= b && b < s, "triangular index");
  return a * s - a * (a - 1) / 2 + (b - a);
}
inline RV symmetric(const Aggregated &ag, const RV &w, size_t s, size_t n) {
  size_t h = s * (s + 1) / 2;
  RV out(4 * h);
  for (int e = 0; e < 4; e++)
    for (size_t a = 0; a < s; a++)
      for (size_t b = a; b < s; b++) {
        auto u = dotR(slice(ag.phi[e], a * n, n), slice(w, b * n, n));
        auto v = dotR(slice(ag.phi[e], b * n, n), slice(w, a * n, n));
        out[e * h + tri(s, a, b)] = scaleR(ra(u, v), iv2);
      }
  for (int e = 0; e < 4; e++) {
    U v = 0;
    for (size_t j = 0; j < s; j++)
      v = add(v, out[e * h + tri(s, j, j)][0]);
    check(v == ag.rhs[e], "honest diagonal aggregated equation");
  }
  return out;
}
inline RV response(const RV &w, const RV &c, size_t n) {
  check(w.size() == n * c.size(), "short response blocks");
  RV z(n);
  for (size_t j = 0; j < c.size(); j++)
    for (size_t a = 0; a < n; a++)
      z[a] = ra(z[a], rm(c[j], w[j * n + a]));
  return z;
}
inline RV combined_phi(const Aggregated &ag, const RV &c, size_t n, int e) {
  RV x(n);
  for (size_t j = 0; j < c.size(); j++)
    for (size_t a = 0; a < n; a++)
      x[a] = ra(x[a], rm(c[j], ag.phi[e][j * n + a]));
  return x;
}
inline Ring principalRHS(const RV &h, const RV &c, int e) {
  size_t s = c.size(), ht = s * (s + 1) / 2;
  Ring out{};
  for (size_t a = 0; a < s; a++)
    for (size_t b = a; b < s; b++) {
      Ring p = rm(c[a], c[b]);
      if (a != b)
        p = scaleR(p, 2);
      out = ra(out, rm(p, h[e * ht + tri(s, a, b)]));
    }
  return out;
}
inline void principalcheck(const Aggregated &ag, const RV &h, const RV &c,
                           const RV &z) {
  for (int e = 0; e < 4; e++)
    check(dotR(combined_phi(ag, c, z.size(), e), z) == principalRHS(h, c, e),
          "terminal principal equation");
}
inline int ellrho(size_t rho) {
  check(rho >= 2, "radix must be at least two");
  int l = 0;
  U p = 1;
  while (p < q) {
    check(p <= UINT64_MAX / rho, "radix multiplication overflow");
    p *= rho;
    l++;
  }
  return l;
}
struct ChildParams {
  Config cfg;
  Key A, B, D;
  RV u1, u2, c;
  Aggregated ag;
  size_t bt, bh;
  int ell;
};
inline System childsystem(const ChildParams &p) {
  size_t n = p.cfg.n, s = p.cfg.s, rho = p.cfg.rho, ht = s * (s + 1) / 2;
  check(n && s && rho >= 2 && p.ell == ellrho(rho), "child configuration");
  check(p.A.rows == 16 && p.A.cols == n && p.B.rows == 6 && p.D.rows == 6 &&
            p.B.cols == p.bt && p.D.cols == p.bh && p.bt == p.ell * 16 * s &&
            p.bh == p.ell * 4 * ht,
        "child key dimensions");
  check(p.c.size() == s && p.u1.size() == 6 && p.u2.size() == 6,
        "child message dimensions");
  for (const auto &phi : p.ag.phi)
    check(phi.size() == s * n, "child coefficient dimensions");
  System sys;
  sys.raw = 64 * (2 * n + p.bt + p.bh);
  append(sys.y, p.u1);
  append(sys.y, p.u2);
  sys.y.resize(sys.y.size() + 64 * (p.A.rows + 4));
  for (auto a : p.ag.rhs)
    sys.y.push_back(a);
  sys.apply = [p, n, s, rho, ht](const RV &x) {
    check(x.size() * 64 >= 64 * (2 * n + p.bt + p.bh), "child witness length");
    auto lo = slice(x, 0, n), hi = slice(x, n, n), bt = slice(x, 2 * n, p.bt),
         bh = slice(x, 2 * n + p.bt, p.bh);
    RV z(n);
    for (size_t j = 0; j < n; j++)
      z[j] = ra(lo[j], scaleR(hi[j], rho));
    RV t = recompose(bt, rho, p.ell), h = recompose(bh, rho, p.ell);
    SV out;
    append(out, p.B.apply(bt));
    append(out, p.D.apply(bh));
    auto az = p.A.apply(z);
    for (size_t j = 0; j < s; j++)
      for (size_t r = 0; r < p.A.rows; r++)
        az[r] = rs(az[r], rm(p.c[j], t[j * p.A.rows + r]));
    append(out, az);
    for (int e = 0; e < 4; e++) {
      auto eq =
          rs(dotR(combined_phi(p.ag, p.c, n, e), z), principalRHS(h, p.c, e));
      out.insert(out.end(), eq.begin(), eq.end());
    }
    for (int e = 0; e < 4; e++) {
      U a = 0;
      for (size_t j = 0; j < s; j++)
        a = add(a, h[e * ht + tri(s, j, j)][0]);
      out.push_back(a);
    }
    return out;
  };
  sys.adj = [p, n, s, rho, ht](const SV &a) {
    size_t pos = 0;
    auto b1 = chunks(a, pos, p.B.rows, true);
    pos += 64 * p.B.rows;
    auto b2 = chunks(a, pos, p.D.rows, true);
    pos += 64 * p.D.rows;
    auto am = chunks(a, pos, p.A.rows, true);
    pos += 64 * p.A.rows;
    auto pe = chunks(a, pos, 4, true);
    pos += 256;
    check(pos + 4 == a.size(), "child transpose rows");
    RV out(2 * n + p.bt + p.bh);
    auto zrow = p.A.transpose(am);
    for (int e = 0; e < 4; e++) {
      auto ph = combined_phi(p.ag, p.c, n, e);
      for (size_t j = 0; j < n; j++)
        zrow[j] = ra(zrow[j], rm(pe[e], ph[j]));
    }
    for (size_t j = 0; j < n; j++) {
      out[j] = zrow[j];
      out[n + j] = scaleR(zrow[j], rho);
    }
    auto btrow = p.B.transpose(b1), bhrow = p.D.transpose(b2);
    for (size_t j = 0; j < s; j++)
      for (size_t r = 0; r < p.A.rows; r++) {
        auto u = rm(am[r], p.c[j]);
        U pow = 1;
        for (int l = 0; l < p.ell; l++) {
          size_t b = (j * p.A.rows + r) * p.ell + l;
          btrow[b] = rs(btrow[b], scaleR(u, pow));
          pow = mul(pow, rho);
        }
      }
    for (int e = 0; e < 4; e++)
      for (size_t i = 0; i < s; i++)
        for (size_t j = i; j < s; j++) {
          Ring cprod = rm(p.c[i], p.c[j]);
          if (i != j)
            cprod = scaleR(cprod, 2);
          Ring u = scaleR(rm(pe[e], cprod), q - 1);
          if (i == j)
            u[0] = add(u[0], a[pos + e]);
          U pow = 1;
          for (int l = 0; l < p.ell; l++) {
            size_t b = (e * ht + tri(s, i, j)) * p.ell + l;
            bhrow[b] = ra(bhrow[b], scaleR(u, pow));
            pow = mul(pow, rho);
          }
        }
    std::copy(btrow.begin(), btrow.end(), out.begin() + 2 * n);
    std::copy(bhrow.begin(), bhrow.end(), out.begin() + 2 * n + p.bt);
    return out;
  };
  return sys;
}
struct TerminalResponse {
  RV tothers, pivot, h, z;
};
inline void verifyterminal(const Config &cfg, const Key &A, const Key &bp,
                           const Aggregated &ag, const TerminalResponse &p,
                           const RV &c) {
  check(cfg.s > 0 && cfg.n > 0 && A.rows == 16 && A.cols == cfg.n &&
            bp.rows == 2 && bp.cols == 256,
        "pivot configuration");
  for (const auto &phi : ag.phi)
    check(phi.size() == cfg.s * cfg.n, "pivot coefficient dimensions");
  check(p.z.size() == cfg.n && energy(p.z) <= cfg.G,
        "terminal response norm/shape");
  check(c.size() == cfg.s && !rzero(c[0]), "terminal challenge");
  check(p.tothers.size() == 16 * (cfg.s - 1) && p.pivot.size() == 2,
        "terminal commitments shape");
  auto t0 = A.apply(p.z);
  for (size_t j = 1; j < cfg.s; j++)
    for (size_t a = 0; a < 16; a++)
      t0[a] = rs(t0[a], rm(c[j], p.tothers[(j - 1) * 16 + a]));
  auto inv = inverseR(c[0]);
  for (auto &v : t0)
    v = rm(inv, v);
  check(bp.apply(digits(t0, 8, 16)) == p.pivot, "pivot binding equation");
  size_t ht = cfg.s * (cfg.s + 1) / 2;
  check(p.h.size() == 4 * ht, "terminal symmetric shape");
  auto h = p.h;
  for (int e = 0; e < 4; e++) {
    check(h[e * ht][0] == 0, "omitted constant must be absent/zero");
    U v = ag.rhs[e];
    for (size_t j = 1; j < cfg.s; j++)
      v = sub(v, h[e * ht + tri(cfg.s, j, j)][0]);
    h[e * ht][0] = v;
  }
  principalcheck(ag, h, c, p.z);
}

} // namespace sonata
