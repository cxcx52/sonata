// Build the square table forward, independently of the public endpoint map.
#include "sonata/epoch.hpp"
#include "fixture.hpp"
#include "sonata/epoch_envelope.hpp"
#include <chrono>
#include <iostream>
#include "epoch_pivot.hpp"
#include "field_replay.hpp"
#include <string>
using namespace sonata;
using namespace sonata::fixture;
namespace ep = sonata::epoch;

static E pw(E x, size_t n) {
  E y = one;
  while (n) {
    if (n & 1)
      y = em(y, x);
    x = sq(x);
    n >>= 1;
  }
  return y;
}
static U penalty(const RV &w, size_t g) {
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

static std::vector<E> forward(const ep::Frozen &f, E rho, E eta, const RV &w,
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

static E fold(std::vector<E> table, const std::vector<E> &r) {
  table.resize(size_t(1) << r.size());
  for (auto t : r) {
    for (size_t i = 0; i < table.size() / 2; ++i)
      table[i] = ea(table[2 * i], em(t, es(table[2 * i + 1], table[2 * i])));
    table.resize(table.size() / 2);
  }
  return table.at(0);
}
template <class F> static void rejects(F run) {
  bool rejected = false;
  try {
    run();
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected, "negative case unexpectedly passed");
}

int main(int argc, char **argv) {
  check(argc == 1 ||
            (argc == 2 && (std::string(argv[1]) == "--target" ||
                           std::string(argv[1]) == "--target-sequential" ||
                           std::string(argv[1]) == "--target-direct")),
        "expected no option, --target, --target-direct or --target-sequential");
  bool target_run = argc == 2;
  auto before = std::chrono::steady_clock::now();
  auto stage = [&](const char *label) {
    if (!target_run)
      return;
    auto now = std::chrono::steady_clock::now();
    std::cout << label << " "
              << std::chrono::duration<double>(now - before).count() << " s"
              << std::endl;
    before = now;
  };
  std::mt19937_64 rng(79);
  auto master = std::make_shared<Master>();
  master->data.resize(target_run ? 44352 : 4 * 256);
  for (auto &r : master->data)
    r = randR(rng);
  Key st{master, 4, 219};
  std::vector<int> arities =
      target_run ? std::vector<int>{125, 125, 124, 124, 124, 124, 124, 124}
                 : std::vector<int>{1, 2};
  ep::Frozen f;
  std::vector<State> objects, sources;
  auto object = [&](const State &s) {
    ep::Object o;
    std::copy_n(s.begin(), 5, o.pub.begin());
    RV v(219);
    for (size_t a = 0; a < 73; ++a)
      for (size_t t = 0; t < 4; ++t)
        for (size_t j = 0; j < 48; ++j)
          put(v, 192 * a + 48 * t + j, (s[a][t] >> j) & 1);
    o.commitment = st.apply(v);
    return o;
  };
  for (size_t i = 0; i < arities.size(); ++i) {
    Fixture inst(arities[i], 73 + 10 * i, false, i ? &objects.back() : nullptr);
    if (i == 0)
      objects.push_back(inst.states.front());
    objects.push_back(inst.states.back());
    ep::Step step;
    step.alpha = {U(inst.k + 17), 3, 4, 5};
    step.gamma = inst.gamma;
    for (int j = 1; j <= inst.k; ++j) {
      sources.push_back(inst.states[j]);
      step.sources.push_back(object(inst.states[j]));
    }
    f.steps.push_back(step);
  }
  for (const auto &s : objects)
    f.states.push_back(object(s));
  objects.insert(objects.end(), sources.begin(), sources.end());
  f.pending.resize(arities.size() + 1, E{9, 1, 0, 0});
  f.aux_commitment.resize(4);
  const auto l = ep::layout(f);
  if (target_run)
    check(l.groups == 66064 && sources.size() == 994,
          "h=8 fixed-capacity profile");
  else
    check(l.groups == 1008 && l.table_size() == 55826, "h=2, unequal arities");
  RV w((l.raw() + 63) / 64);
  size_t g = 0;
  for (size_t i = 0; i < objects.size(); ++i) {
    bool fresh = i >= f.states.size();
    for (size_t a = 5; a < (fresh ? 69u : 73u); ++a)
      for (size_t t = 0; t < (fresh ? 1u : 4u); ++t) {
        U v = objects[i][a][t];
        size_t p = 55 * g++;
        for (size_t j = 0; j < 48; ++j)
          put(w, p + j, (v >> j) & 1);
        int z = 42 - __builtin_popcountll(v >> 6), e = z == 0;
        int mask = helper_mask((z - 1 + e) / 2);
        for (int j = 0; j < 5; ++j)
          put(w, p + 48 + j, (mask >> j) & 1);
        put(w, p + 53, e);
        put(w, p + 54, e * ((v >> 2) & 1));
      }
  }
  put(w, l.raw() - 1, 1);
  if (argc == 2 && std::string(argv[1]) == "--target-sequential") {
    // Reuse the exact same public state/source commitments and canonical
    // encodings; each standalone settlement needs its own auxiliary tag.
    stage("shared fixture, canonical witness and state/source commitments");
    size_t next_source = f.states.size(), round_sum = 0, table_sum = 0,
           row_sum = 0;
    std::array<double, 5> seconds{};
    auto measure = [&](size_t k) {
      auto now = std::chrono::steady_clock::now();
      seconds[k] += std::chrono::duration<double>(now - before).count();
      before = now;
    };
    for (size_t i = 0; i < f.steps.size(); ++i) {
      ep::Frozen single{{f.states[i], f.states[i + 1]},
                        {f.pending[i], f.pending[i + 1]},
                        {f.steps[i]},
                        RV(4)};
      auto local = ep::layout(single);
      std::vector<size_t> selected{i, i + 1};
      for (size_t j = 0; j < arities[i]; ++j)
        selected.push_back(next_source++);
      RV local_w((local.raw() + 63) / 64);
      for (size_t j = 0; j < selected.size(); ++j) {
        size_t src = l.blocks[selected[j]].first * 55,
               dest = 55 * local.blocks[j].first;
        for (size_t z = 0; z < 55 * (j < 2 ? 272u : 64u); ++z)
          put(local_w, dest + z, cv(w, src + z));
      }
      put(local_w, local.raw() - 1, 1);
      Key small_aux{master, 4, (7 * local.groups + 63) / 64};
      RV av(small_aux.cols);
      for (size_t j = 0; j < 7 * local.groups; ++j)
        put(av, j, cv(local_w, 55 * (j / 7) + 48 + j % 7));
      single.aux_commitment = small_aux.apply(av);
      measure(0);
      E rr = randomE(rng), ee = randomE(rng);
      auto table = forward(single, rr, ee, local_w, true);
      measure(1);
      auto tr = prove(std::move(table), E{mul(local.groups, KAPPA), 0, 0, 0},
                      13 + i, true);
      measure(2);
      auto sys =
          ep::initialsystem(single, rr, ee, tr, master, small_aux.cols, true);
      measure(3);
      check(sys.apply(local_w) == sys.y, "standalone direct-settlement root");
      measure(4);
      round_sum += tr.ch.size();
      row_sum += sys.y.size();
      table_sum += local.table_size();
    }
    constexpr const char *labels[] = {
        "split witness and per-step auxiliary commitments",
        "dense forward table and independent equation sum", "square sumcheck",
        "public root and endpoint operator", "active root apply"};
    for (size_t j = 0; j < seconds.size(); ++j)
      std::cout << labels[j] << " " << seconds[j] << " s\n";
    std::cout << "PASS eight h=1 direct settlements n=994 T_sum=" << table_sum
              << " active_rows_sum=" << row_sum << " rounds_sum=" << round_sum
              << "\n";
    return 0;
  }
  Key aux{master, 4, target_run ? size_t(7228) : size_t(128)};
  RV av(aux.cols);
  for (size_t j = 0; j < 7 * g; ++j)
    put(av, j, cv(w, 55 * (j / 7) + 48 + j % 7));
  f.aux_commitment = aux.apply(av);
  stage("fixture and original commitments");
  RV arbitrary;
  if (!target_run) {
    arbitrary = w;
    for (auto &r : arbitrary)
      r = randR(rng);
  }
  E rho = randomE(rng), eta = randomE(rng);
  for (bool direct : {false, true}) {
    if (target_run && direct != (std::string(argv[1]) == "--target-direct"))
      continue;
    auto table = forward(f, rho, eta, w, direct);
    stage("dense forward table and independent equation sum");
    E claim{mul(g, KAPPA), 0, 0, 0};
    auto tr = prove(std::move(table), claim, 13, true);
    stage("square sumcheck");
    auto sys = ep::initialsystem(f, rho, eta, tr, master, aux.cols, direct);
    stage("public root and endpoint operator");
    if (target_run) {
      constexpr Config cfg{21,      2705,      64,   3634401, 3635520,
                           3105761, 795074815, 1174, 22};
      const size_t active_rows = sys.y.size(), active_raw = sys.raw;
      auto physical = align(ep::embed_prefix(std::move(sys), cfg, 264036), cfg);
      w.resize(cfg.N / 64);
      put(w, active_raw, 17);
      check(physical.apply(w) == physical.y,
            "aligned full-capacity root, with ignored-tail perturbation");
      stage("aligned root apply including tail perturbation");
      std::cout << "PASS h=8 n=994 V=" << l.groups << " T=" << l.table_size()
                << " raw=" << active_raw << " active_rows=" << active_rows
                << " aligned_rows=" << physical.y.size()
                << " rounds=" << tr.ch.size()
                << " drop_pending_equations=" << direct << '\n';
      return 0;
    }
    check(sys.apply(w) == sys.y, "frozen-epoch honest public root");
    if (direct) check_frozen_honest_handoff(f,rho,eta,tr,w,master);
    else check_epoch_pivot(sys,w,master);
    checkadj(sys, arbitrary, rng);
    auto ar = forward(f, rho, eta, arbitrary, direct);
    E target = es(fold(ar, tr.ch), tr.terminal);
    SV rows = sys.apply(arbitrary);
    for (size_t i = 0; i < 4; ++i) {
      size_t pos = rows.size() - 5 + i;
      check(sub(rows[pos], sys.y[pos]) == target[i],
            "nonbit epoch original terminal residual");
    }
    for (size_t idx : {size_t(0), l.blocks[1].first * 55,
                       l.blocks.back().first * 55, size_t(49), l.raw() - 1}) {
      RV bad = w;
      put(bad, idx, add(cv(bad, idx), 1));
      check(sys.apply(bad) != sys.y, "state/source/aux or fixed-one mutation");
    }
    auto bad = f;
    bad.states[1].pub[0][0] = add(bad.states[1].pub[0][0], 1);
    rejects([&] {
      ep::initialsystem(bad, rho, eta, tr, master, aux.cols, direct);
    });
    bad = f;
    bad.pending.back()[0] = add(bad.pending.back()[0], 1);
    if (direct)
      rejects([&] {
        ep::initialsystem(bad, rho, eta, tr, master, aux.cols, true);
      });
    else {
      auto altered = ep::initialsystem(bad, rho, eta, tr, master, aux.cols);
      check(altered.apply(w) != altered.y, "post-freeze pending mutation");
      auto mismatch = forward(bad, rho, eta, arbitrary, false);
      check(fold(mismatch, tr.ch) != fold(ar, tr.ch),
            "pending contributes to independent endpoint");
    }
  }
  std::cout << "PASS two-step epoch, arities 1/2, both equation schedules, "
               "arbitrary field inputs and original sumcheck endpoint\n";
}
