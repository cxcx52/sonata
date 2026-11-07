#pragma once
#include "ring.hpp"
namespace sonata {
struct PublicFront {
  int k;
  E alpha, rho, gamma, eta, prev_pending, pending;
  std::vector<std::array<E, 5>> pub;
  std::vector<RV> cs;
  RV caux;
  Transcript field;
};
struct AffineFront {
  std::vector<E> row;
  E off{};
  std::vector<std::array<int, 4>> map;
};
inline AffineFront fieldadj(const PublicFront &f, const std::vector<E> &ew) {
  size_t V = 64 * f.k + 544;
  AffineFront o;
  o.row.resize(55 * V);
  o.map.resize((f.k + 2) * 73);
  for (auto &a : o.map)
    a.fill(-1);
  size_t gs = 0;
  for (int s = 0; s < f.k + 2; s++)
    for (int a = 5; a < 73; a++) {
      bool fresh = s > 0 && s <= f.k;
      if (fresh && a >= 69)
        continue;
      for (int t = 0; t < (fresh ? 1 : 4); t++)
        o.map[73 * s + a][t] = gs++;
    }
  check(gs == V, "public layout V");
  size_t pos = 0;
  for (size_t g = 0; g < V; g++) {
    E prefix{};
    for (int i = 0; i < 55; i++) {
      E rt{SQRT_D[i][0], SQRT_D[i][1], SQRT_D[i][2], SQRT_D[i][3]};
      E a = em(ew.at(pos++), rt);
      o.off = ea(o.off, sc(a, BETA[i]));
      if (i >= 6 && i < 48) {
        o.row[55 * g + i] = ea(o.row[55 * g + i], ea(a, prefix));
        prefix = ea(prefix, sc(a, high_inv[i - 6]));
      } else
        for (int j = i; j < 55; j++)
          if (LDL_L[j][i])
            o.row[55 * g + j] = ea(o.row[55 * g + j], sc(a, LDL_L[j][i]));
    }
    for (int j = 0; j < 5; j++)
      o.row[55 * g + 48 + j] =
          ea(o.row[55 * g + 48 + j], sc(prefix, 2 * HWEIGHTS[j]));
    o.row[55 * g + 53] = es(o.row[55 * g + 53], prefix);
  }
  auto decadj = [&](int s, int c, E alpha) {
    for (int t = 0; t < 4; t++) {
      int g = o.map[73 * s + c][t];
      E basis{};
      basis[t] = 1;
      E a = em(alpha, basis);
      if (g < 0) {
        U z = c < 5 ? f.pub[s][c][t] : 0;
        o.off = ea(o.off, sc(a, z));
        continue;
      }
      for (int j = 0; j < 48; j++) {
        o.row[55 * g + j] = ea(o.row[55 * g + j], a);
        a = sc(a, 2);
      }
    }
  };
  std::vector<Lin> linears;
  E fpower = one, spower = one, wt = one, lconst{};
  auto next = [&](bool fresh) {
    if (fresh) {
      fpower = em(fpower, f.gamma);
      return em(f.eta, fpower);
    }
    spower = em(spower, f.rho);
    return em(f.eta, spower);
  };
  for (int s = 0; s < f.k + 2; s++)
    for (int j = 0; j < 4; j++) {
      bool fresh = s > 0 && s <= f.k;
      wt = next(fresh);
      for (int l = 0; l < 8; l++) {
        E u = ew.at(pos++), v = sc(ew.at(pos++), IOTA);
        decadj(s, 5 + 16 * j + 2 * l, sc(em(wt, ea(u, v)), iv2));
        decadj(s, 6 + 16 * j + 2 * l, sc(es(u, v), iv2));
      }
      lconst = es(lconst, em(wt, em(f.pub[s][0], f.pub[s][1 + j])));
      if (!fresh)
        linears.push_back({es(zero, wt), s, 69 + j});
    }
  std::vector<E> powers(f.k);
  E p = one;
  for (auto &v : powers) {
    p = em(p, f.alpha);
    v = p;
  }
  for (int a = 5; a < 69; a++) {
    wt = next(false);
    linears.push_back({wt, f.k + 1, a});
    linears.push_back({es(zero, wt), 0, a});
    for (int s = 1; s <= f.k; s++)
      linears.push_back({es(zero, em(wt, powers[s - 1])), s, a});
  }
  E u = ew.at(pos++), v = sc(ew.at(pos++), IOTA), alpha = sc(ea(u, v), iv2);
  o.off = ea(o.off, sc(es(u, v), iv2));
  o.off = ea(o.off, em(alpha, lconst));
  for (auto &r : linears)
    decadj(r.s, r.c, em(alpha, r.w));
  return o;
}
inline void frontvalidate(const PublicFront &f) {
  check(f.k >= 1 && f.k <= 1024, "unsupported root arity");
  check(f.caux.size() == 4 &&
            std::all_of(f.cs.begin(), f.cs.end(),
                        [](const RV &v) { return v.size() == 4; }),
        "commitment shape");
  check(f.pub.size() == size_t(f.k + 2) && f.cs.size() == f.pub.size(),
        "front shapes");
  for (int s = 1; s <= f.k; s++) {
    check(f.pub[s][0] == one, "fresh u");
    for (int a = 1; a < 5; a++)
      for (int t = 1; t < 4; t++)
        check(f.pub[s][a][t] == 0, "fresh base field");
  }
  for (int a = 0; a < 5; a++) {
    E x = f.pub[0][a], p = one;
    for (int s = 1; s <= f.k; s++) {
      p = em(p, f.alpha);
      x = ea(x, em(p, f.pub[s][a]));
    }
    check(x == f.pub[f.k + 1][a], "folded public coordinates");
  }
  check(f.gamma != zero && f.rho != zero && f.eta != zero,
        "nonzero batch challenges");
  size_t table = 55 * (64 * f.k + 544) + 64 * (f.k + 2) + 2, P = 1, ell = 0;
  while (P < table)
    P *= 2, ell++;
  check(f.field.ch.size() == ell, "sumcheck rounds");
  E claim = ea(E{mul(64 * f.k + 544, KAPPA), 0, 0, 0},
               em(f.eta, es(f.pending, f.prev_pending)));
  check(verify(f.field, claim, f.field.terminal),
        "field transcript recurrence");
}
inline System initialsystem(const PublicFront &f,
                            std::shared_ptr<const Master> master) {
  frontvalidate(f);
  auto ew = eqweights(f.field.ch);
  auto af = std::make_shared<AffineFront>(fieldadj(f, ew));
  size_t V = 64 * f.k + 544, raw = 55 * V + 1, auxcols = (7 * V + 63) / 64;
  Key st{master, 4, 219}, aux{master, 4, auxcols};
  std::vector<std::vector<std::pair<size_t, size_t>>> ins(f.k + 3);
  std::vector<RV> fixed(f.k + 2, RV(219));
  for (int s = 0; s < f.k + 2; s++)
    for (int a = 0; a < 73; a++)
      for (int t = 0; t < 4; t++)
        for (int j = 0; j < 48; j++) {
          size_t phys = 192 * a + 48 * t + j;
          int g = af->map[73 * s + a][t];
          if (g >= 0)
            ins[s].push_back({phys, 55 * size_t(g) + j});
          else if (a < 5)
            put(fixed[s], phys, (f.pub[s][a][t] >> j) & 1);
        }
  for (size_t g = 0; g < V; g++)
    for (int t = 0; t < 7; t++)
      ins[f.k + 2].push_back({7 * g + t, 55 * g + 48 + t});
  System sys;
  sys.raw = raw;
  for (int s = 0; s < f.k + 2; s++) {
    RV fix = st.apply(fixed[s]), active(4);
    for (int a = 0; a < 4; a++)
      active[a] = rs(f.cs[s][a], fix[a]);
    append(sys.y, active);
  }
  append(sys.y, f.caux);
  E y = es(f.field.terminal, af->off);
  for (U a : y)
    sys.y.push_back(a);
  sys.y.push_back(1);
  sys.apply = [ins, af, st, aux, raw, V, k = f.k](const RV &w) {
    SV y;
    for (int s = 0; s < k + 3; s++) {
      Key key = s == k + 2 ? aux : st;
      RV v(key.cols);
      for (auto [p, a] : ins[s])
        put(v, p, cv(w, a));
      append(y, key.apply(v));
    }
    E t{};
    for (size_t i = 0; i < 55 * V; i++)
      t = ea(t, sc(af->row[i], cv(w, i)));
    for (U a : t)
      y.push_back(a);
    y.push_back(cv(w, raw - 1));
    return y;
  };
  sys.adj = [ins, af, st, aux, raw, V, k = f.k](const SV &a) {
    RV coeff((raw + 63) / 64);
    size_t p = 0;
    for (int s = 0; s < k + 3; s++) {
      Key key = s == k + 2 ? aux : st;
      auto v = key.transpose(chunks(a, p, 4, true));
      p += 256;
      for (auto [phys, act] : ins[s]) {
        auto d = dual(v[phys / 64]);
        put(coeff, act, add(cv(coeff, act), d[phys % 64]));
      }
    }
    for (size_t j = 0; j < 55 * V; j++) {
      U z = 0;
      for (int t = 0; t < 4; t++)
        z = add(z, mul(a[p + t], af->row[j][t]));
      put(coeff, j, add(cv(coeff, j), z));
    }
    put(coeff, raw - 1, add(cv(coeff, raw - 1), a[p + 4]));
    for (auto &r : coeff)
      r = dual(r);
    return coeff;
  };
  return sys;
}

} // namespace sonata
