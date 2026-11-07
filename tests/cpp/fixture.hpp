#pragma once
#include "sonata/frontend.hpp"
#include <random>
namespace sonata::fixture {
static E randomE(std::mt19937_64 &r) {
  E a;
  for (auto &v : a)
    v = r() % q;
  return a;
}
static int helper_mask(int t) {
  int best = -1;
  for (int mask = 0; mask < 32; mask++) {
    int v = 0;
    for (int j = 0; j < 5; j++)
      if ((mask >> j) & 1)
        v += HWEIGHTS[j];
    if (v == t && (best < 0 || __builtin_popcount(unsigned(mask)) <
                                   __builtin_popcount(unsigned(best))))
      best = mask;
  }
  check(best >= 0, "aux count representability");
  return best;
}
using State = std::array<E, 73>;
struct Product {
  E w;
  int sa, ca, sb, cb;
};
struct Fixture {
  int k;
  std::vector<State> states;
  std::vector<E> weights;
  std::vector<std::pair<int, int>> coords;
  std::vector<int> components;
  std::vector<std::array<int, 4>> gmap;
  std::vector<U> scalarvalues;
  std::vector<uint8_t> bits;
  std::vector<Product> products;
  std::vector<Lin> linears;
  E lconst{}, rho{}, gamma{}, eta{}, pending{}, linear{};
  size_t constraints = 0;
  int invalid_sources = 0;
  Fixture(int kk, U seed, bool invalid = false,
          const State *old_state = nullptr)
      : k(kk), states(kk + 2), gmap((kk + 2) * 73) {
    for (auto &a : gmap)
      a.fill(-1);
    std::mt19937_64 rng(seed);
    for (int s = 0; s <= k; s++) {
      auto &v = states[s];
      v[0] = s ? one : randomE(rng);
      for (int j = 0; j < 4; j++) {
        E d{};
        for (int l = 0; l < 8; l++) {
          int a = 5 + 16 * j + 2 * l;
          v[a] = s ? E{rng() % q, 0, 0, 0} : randomE(rng);
          v[a + 1] = s ? E{rng() % q, 0, 0, 0} : randomE(rng);
          d = ea(d, em(v[a], v[a + 1]));
        }
        v[1 + j] = s ? d : randomE(rng);
        v[69 + j] = es(d, em(v[0], v[1 + j]));
      }
    }
    // Reuse the previous node's exact output state, rather than another random
    // state.
    if (old_state)
      states[0] = *old_state;
    // Make one fresh instance invalid while retaining a valid fold/output
    // relation. In the commitment game this public source input is fixed before
    // gamma.
    if (invalid)
      states[1][1][0] = add(states[1][1][0], 1);
    E alpha{U(k + 17), 3, 4, 5};
    weights.resize(k);
    E power = one;
    for (int i = 0; i < k; i++) {
      power = em(power, alpha);
      weights[i] = power;
    }
    for (int a = 0; a < 69; a++) {
      states[k + 1][a] = states[0][a];
      for (int i = 0; i < k; i++)
        states[k + 1][a] =
            ea(states[k + 1][a], em(weights[i], states[i + 1][a]));
    }
    for (int j = 0; j < 4; j++) {
      E d{};
      for (int l = 0; l < 8; l++)
        d = ea(d, em(states[k + 1][5 + 16 * j + 2 * l],
                     states[k + 1][6 + 16 * j + 2 * l]));
      states[k + 1][69 + j] = es(d, em(states[k + 1][0], states[k + 1][1 + j]));
    }
    for (int s = 0; s < k + 2; s++)
      for (int a = 5; a < 73; a++) {
        bool fresh = s > 0 && s <= k;
        if (fresh && a >= 69)
          continue;
        for (int t = 0; t < (fresh ? 1 : 4); t++) {
          int g = scalarvalues.size();
          gmap[73 * s + a][t] = g;
          scalarvalues.push_back(states[s][a][t]);
          coords.push_back({s, a});
          components.push_back(t);
        }
      }
    check(scalarvalues.size() == size_t(64 * k + 544),
          "public-specialized coefficient count");
    bits.resize(55 * scalarvalues.size());
    for (size_t g = 0; g < scalarvalues.size(); g++) {
      U x = scalarvalues[g];
      auto b = bits.data() + 55 * g;
      for (int j = 0; j < 48; j++)
        b[j] = (x >> j) & 1;
      int z = 42 - __builtin_popcountll(x >> 6), e = z == 0,
          t = (z - 1 + e) / 2;
      int mask = helper_mask(t);
      for (int j = 0; j < 5; j++)
        b[48 + j] = (mask >> j) & 1;
      b[53] = e;
      b[54] = e * b[2];
    }
    // Challenge chronology: source commitments -> gamma -> pending claim ->
    // rho, eta. A real adversarial-security proof must bind source openings
    // before gamma.
    gamma = randomE(rng);
    if (gamma == zero)
      gamma = one;
    rho = randomE(rng);
    if (rho == zero)
      rho = one;
    eta = randomE(rng);
    if (eta == zero)
      eta = one;
    E fresh_power = one, strict_power = one;
    auto next = [&](bool fresh) {
      constraints++;
      if (fresh) {
        fresh_power = em(fresh_power, gamma);
        return em(eta, fresh_power);
      }
      strict_power = em(strict_power, rho);
      return em(eta, strict_power);
    };
    for (int s = 0; s < k + 2; s++)
      for (int j = 0; j < 4; j++) {
        bool fresh = (s > 0 && s <= k);
        E w = next(fresh), res{};
        for (int l = 0; l < 8; l++) {
          int a = 5 + 16 * j + 2 * l;
          products.push_back({w, s, a, s, a + 1});
          res = ea(res, em(states[s][a], states[s][a + 1]));
        }
        E pub = em(states[s][0], states[s][1 + j]);
        lconst = es(lconst, em(w, pub));
        res = es(res, pub);
        if (!fresh) {
          linears.push_back({es(zero, w), s, 69 + j});
          res = es(res, states[s][69 + j]);
        }
        if (fresh) {
          pending = ea(pending, em(fresh_power, res));
          if (res != zero)
            invalid_sources++;
        } else
          check(res == zero, "strict old/output quadratic relation");
      }
    for (int a = 5; a < 69; a++) {
      E w = next(false);
      linears.push_back({w, k + 1, a});
      linears.push_back({es(zero, w), 0, a});
      E res = es(states[k + 1][a], states[0][a]);
      for (int s = 1; s <= k; s++) {
        E v = em(w, weights[s - 1]);
        linears.push_back({es(zero, v), s, a});
        res = es(res, em(weights[s - 1], states[s][a]));
      }
      check(res == zero, "private fold equation");
    }
    linear = lconst;
    for (auto &r : linears)
      linear = ea(linear, em(r.w, states[r.s][r.c]));
    check(products.size() == size_t(32 * (k + 2)), "nonlinear products");
    check(constraints == size_t(4 * k + 72), "semantic constraints");
  }
  std::vector<E> table() const {
    std::vector<E> t;
    t.reserve(55 * scalarvalues.size() + 2 * products.size() + 2);
    for (size_t g = 0; g < scalarvalues.size(); g++) {
      auto b = bits.data() + 55 * g;
      U suffix = 0;
      for (int j = 6; j < 48; j++)
        suffix += b[j];
      for (int j = 0; j < 5; j++)
        suffix += U(b[48 + j]) * 2 * HWEIGHTS[j];
      suffix = sub(suffix, b[53]);
      for (int i = 0; i < 55; i++) {
        U v = BETA[i];
        if (i >= 6 && i < 48) {
          suffix = sub(suffix, b[i]);
          v = add(v, add(b[i], mul(high_inv[i - 6], suffix)));
        } else
          for (int j = i; j < 55; j++)
            if (b[j])
              v = add(v, LDL_L[j][i]);
        E rt{SQRT_D[i][0], SQRT_D[i][1], SQRT_D[i][2], SQRT_D[i][3]};
        t.push_back(sc(rt, v));
      }
    }
    U half = (q + 1) / 2;
    for (auto &p : products) {
      E wa = em(p.w, states[p.sa][p.ca]), b = states[p.sb][p.cb];
      t.push_back(sc(ea(wa, b), half));
      t.push_back(sc(es(wa, b), mul(IOTA, half)));
    }
    t.push_back(sc(ea(linear, one), half));
    t.push_back(sc(es(linear, one), mul(IOTA, half)));
    return t;
  }
};
static Transcript prove(std::vector<E> a, E claim, U seed,
                        bool compact = false) {
  Transcript tr;
  std::mt19937_64 rng(seed);
  size_t P = 1;
  while (P < a.size())
    P *= 2;
  if (!compact)
    a.resize(P);
  while (a.size() > 1) {
    size_t n = (a.size() + 1) / 2;
    E g0{}, g2{};
    for (size_t i = 0; i < n; i++) {
      E right = 2 * i + 1 < a.size() ? a[2 * i + 1] : zero;
      E d = es(right, a[2 * i]);
      g0 = ea(g0, sq(a[2 * i]));
      g2 = ea(g2, sq(d));
      tr.squares += 2;
    }
    E g1 = es(es(claim, sc(g0, 2)), g2);
    E r = randomE(rng); // test coins; sampled only after round coefficients
    tr.g0.push_back(g0);
    tr.g2.push_back(g2);
    tr.ch.push_back(r);
    claim = ea(g0, em(r, ea(g1, em(r, g2))));
    for (size_t i = 0; i < n; i++) {
      E right = 2 * i + 1 < a.size() ? a[2 * i + 1] : zero;
      a[i] = ea(a[2 * i], em(r, es(right, a[2 * i])));
      tr.foldmults++;
    }
    a.resize(n);
  }
  tr.terminal = a[0];
  check(sq(tr.terminal) == claim, "prover's terminal claim");
  tr.mults = tr.squares + tr.foldmults;
  return tr;
}
static U randomQ(std::mt19937_64 &r) {
  U x;
  do {
    x = r() & mask;
  } while (x >= q);
  return x;
}
static Ring randR(std::mt19937_64 &r) {
  Ring a;
  for (auto &x : a)
    x = randomQ(r);
  return a;
}
static void checkadj(const System &sys, const RV &w, std::mt19937_64 &r) {
  SV a(sys.y.size());
  for (auto &x : a)
    x = randomQ(r);
  auto y = sys.apply(w);
  U d = 0;
  for (size_t i = 0; i < y.size(); i++)
    d = add(d, mul(a[i], y[i]));
  check(ctDot(sys.adj(a), w) == d, "scalar/ring adjoint identity");
}

} // namespace sonata::fixture
