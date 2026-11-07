#pragma once
#include "fixture.hpp"
#include "sonata/epoch.hpp"

namespace sonata::fixture {
namespace ep = sonata::epoch;

// Independent prover-side square table used by the small and full-profile fixtures.
inline E pw(E x, size_t n) {
  E y = one;
  while (n) {
    if (n & 1)
      y = em(y, x);
    x = sq(x);
    n >>= 1;
  }
  return y;
}
inline U penalty(const RV &w, size_t g) {
  auto x = [&](size_t j) { return cv(w, 55 * g + j); };
  U binary = 0, z = 41;
  for (size_t j = 0; j < 55; ++j)
    binary = add(binary, mul(x(j), sub(x(j), 1)));
  for (size_t j = 6; j < 48; ++j)
    z = sub(z, x(j));
  for (size_t j = 0; j < 5; ++j)
    z = sub(z, mul(2 * HWEIGHTS[j], x(48 + j)));
  U e = x(53), m = x(54), b = x(2);
  z = add(z, e);
  U f = add(mul(e, b), mul(3, m));
  f = sub(f, mul(2, mul(m, add(e, b))));
  f = add(f, mul(e, add(x(3), add(x(4), x(5)))));
  f = add(f, mul(m, add(x(0), x(1))));
  return add(mul(8, binary), add(mul(z, sub(z, 1)), f));
}

inline std::vector<E> forward(const ep::Frozen &f, E rho, E eta, const RV &w,
                              bool direct) {
  const size_t h = f.steps.size();
  std::vector<State> objs;
  std::vector<E> table;
  table.reserve(ep::layout(f).table_size());
  size_t group = 0;
  U energy_sum = 0;
  auto decode = [&](const ep::Object &o, bool fresh) {
    State s{};
    std::copy(o.pub.begin(), o.pub.end(), s.begin());
    for (size_t slot = 5; slot < (fresh ? 69u : 73u); ++slot)
      for (size_t t = 0; t < (fresh ? 1u : 4u); ++t) {
        for (int j = 47; j >= 0; --j)
          s[slot][t] = add(mul(s[slot][t], 2), cv(w, 55 * group + j));
        for (size_t i = 0; i < 55; ++i) {
          U form = BETA[i];
          for (size_t j = 0; j < 55; ++j)
            form = add(form, mul(LDL_L[j][i], cv(w, 55 * group + j)));
          table.push_back(sc(
              {SQRT_D[i][0], SQRT_D[i][1], SQRT_D[i][2], SQRT_D[i][3]}, form));
        }
        energy_sum = add(energy_sum, penalty(w, group++));
      }
    objs.push_back(s);
  };
  for (const auto &o : f.states)
    decode(o, false);
  for (const auto &t : f.steps)
    for (const auto &o : t.sources)
      decode(o, true);
  std::vector<std::array<E, 4>> weights(objs.size());
  std::vector<E> equations;
  E linear{};
  auto residual = [&](size_t obj, size_t j) {
    const State &s = objs[obj];
    E value = es(zero, em(s[0], s[1 + j]));
    for (size_t z = 0; z < 8; ++z)
      value = ea(value, em(s[5 + 16 * j + 2 * z], s[6 + 16 * j + 2 * z]));
    if (obj <= h)
      value = es(value, s[69 + j]);
    return value;
  };
  for (size_t obj = 0; obj <= h; ++obj)
    for (size_t j = 0; j < 4; ++j) {
      equations.push_back(residual(obj, j));
      weights[obj][j] = em(eta, pw(rho, equations.size()));
    }
  size_t first = h + 1;
  for (size_t i = 0; i < h; ++i) {
    const auto &step = f.steps[i];
    size_t k = step.sources.size();
    size_t start = equations.size();
    E pending = es(f.pending[i + 1], f.pending[i]);
    for (size_t j = 0; j < 4 * k; ++j) {
      E res = residual(first + j / 4, j % 4);
      equations.push_back(res);
      pending = es(pending, em(pw(step.gamma, j + 1), res));
      E combined = pw(rho, start + j + 1);
      if (!direct)
        combined = es(combined,
                      em(pw(rho, start + 4 * k + 65), pw(step.gamma, j + 1)));
      weights[first + j / 4][j % 4] = em(eta, combined);
    }
    for (size_t a = 5; a < 69; ++a) {
      E fold = es(objs[i + 1][a], objs[i][a]);
      for (size_t j = 0; j < k; ++j)
        fold = es(fold, em(pw(step.alpha, j + 1), objs[first + j][a]));
      equations.push_back(fold);
      linear = ea(linear, em(em(eta, pw(rho, equations.size())), fold));
    }
    if (!direct) {
      equations.push_back(pending);
      linear = ea(linear, em(em(eta, pw(rho, equations.size())),
                             es(f.pending[i + 1], f.pending[i])));
    }
    first += k;
  }
  auto pair = [&](E a, E b) {
    table.push_back(sc(ea(a, b), iv2));
    table.push_back(sc(es(a, b), mul(iv2, IOTA)));
  };
  for (size_t obj = 0; obj < objs.size(); ++obj)
    for (size_t j = 0; j < 4; ++j) {
      const auto &s = objs[obj];
      E wt = weights[obj][j];
      for (size_t z = 0; z < 8; ++z)
        pair(em(wt, s[5 + 16 * j + 2 * z]), s[6 + 16 * j + 2 * z]);
      E aff = em(s[0], s[1 + j]);
      if (obj <= h)
        aff = ea(aff, s[69 + j]);
      linear = es(linear, em(wt, aff));
    }
  pair(linear, one);
  E square_sum{}, expected{add(mul(group, KAPPA), energy_sum), 0, 0, 0};
  for (const auto &v : table)
    square_sum = ea(square_sum, sq(v));
  for (size_t j = 0; j < equations.size(); ++j)
    expected = ea(expected, em(eta, em(pw(rho, j + 1), equations[j])));
  check(
      square_sum == expected,
      "square table vs integer encoding polynomial and all semantic equations");
  return table;
}

struct WindowOpening {
  ep::Frozen log;
  RV witness;
};

inline std::vector<int> target_window_arities() {
  return {125, 125, 124, 124, 124, 124, 124, 124};
}

inline ep::Object commit_state(const Key& st, const State& s) {
  ep::Object o;
  std::copy_n(s.begin(), 5, o.pub.begin());
  RV v(219);
  for (size_t a = 0; a < 73; ++a)
    for (size_t t = 0; t < 4; ++t)
      for (size_t j = 0; j < 48; ++j)
        put(v, 192 * a + 48 * t + j, (s[a][t] >> j) & 1);
  o.commitment = st.apply(v);
  return o;
}

inline RV window_witness(const ep::Frozen& f, const std::vector<State>& objects) {
  const auto l = ep::layout(f);
  RV w((l.raw() + 63) / 64);
  size_t g = 0;
  for (size_t i = 0; i < objects.size(); ++i) {
    bool fresh = i >= f.states.size();
    for (size_t a = 5; a < (fresh ? 69u : 73u); ++a)
      for (size_t t = 0; t < (fresh ? 1u : 4u); ++t) {
        const auto group = canonical_group(objects[i][a][t]);
        size_t p = 55 * g++;
        for (size_t j = 0; j < group.size(); ++j)
          put(w, p + j, group[j]);
      }
  }
  check(g == l.groups, "canonical group count");
  put(w, l.raw() - 1, 1);
  return w;
}

// The history and its public anchor are reproducible inputs to certification.
inline WindowOpening make_window_opening(std::shared_ptr<const Master> master,
                                         const std::vector<int>& arities,
                                         E initial_pending) {
  Key st{master, 4, 219};
  ep::Frozen f;
  f.pending.push_back(initial_pending);
  std::vector<State> objects, sources;
  for (size_t i = 0; i < arities.size(); ++i) {
    Fixture inst(arities[i], 73 + 10 * i, false, i ? &objects.back() : nullptr);
    if (i == 0)
      objects.push_back(inst.states.front());
    objects.push_back(inst.states.back());
    ep::Step step;
    step.alpha = inst.weights.front();
    step.gamma = inst.gamma;
    for (int j = 1; j <= inst.k; ++j) {
      sources.push_back(inst.states[j]);
      step.sources.push_back(commit_state(st, inst.states[j]));
    }
    check(inst.invalid_sources == 0 && inst.pending == zero, "expected valid source data");
    f.pending.push_back(ea(f.pending.back(), inst.pending));
    f.steps.push_back(std::move(step));
  }
  for (const auto &s : objects)
    f.states.push_back(commit_state(st, s));
  objects.insert(objects.end(), sources.begin(), sources.end());
  f.aux_commitment.resize(4);
  auto w = window_witness(f, objects);
  return {std::move(f), std::move(w)};
}

inline void commit_window_aux(ep::Frozen& log, const RV& witness,
                              const Key& aux, size_t groups) {
  RV opening(aux.cols);
  for (size_t j = 0; j < 7 * groups; ++j)
    put(opening, j, cv(witness, 55 * (j / 7) + 48 + j % 7));
  log.aux_commitment = aux.apply(opening);
}

} // namespace sonata::fixture
