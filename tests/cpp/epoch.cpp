// Build the square table forward, independently of the public endpoint map.
#include "epoch_fixture.hpp"
#include "sonata/epoch_envelope.hpp"
#include <chrono>
#include <iostream>
#include "epoch_pivot.hpp"
#include "field_replay.hpp"
#include <string>
using namespace sonata;
using namespace sonata::fixture;
namespace ep = sonata::epoch;

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
  const auto arities = target_run ? target_window_arities() : std::vector<int>{1, 2};
  auto [f, w] = make_window_opening(master, arities, E{9, 1, 0, 0});
  const auto l = ep::layout(f);
  const size_t g = l.groups;
  if (target_run)
    check(g == 66064 && l.blocks.size() - f.states.size() == 994, "h=8 fixed-capacity profile");
  else
    check(g == 1008 && l.table_size() == 55826, "h=2, unequal arities");
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
      commit_window_aux(single, local_w, small_aux, local.groups);
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
  commit_window_aux(f, w, aux, g);
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
