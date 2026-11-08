#pragma once

#include "sonata/backend.hpp"
#include <chrono>
#include <iostream>
#include <string_view>

namespace sonata::example {

// Rows of parameters/ledger-two-point.json, with the radices from compute.py.
// For this profile the response energy cap is 256*S-1, or 256*S at the pivot.
inline constexpr std::array<Config, 5> layers{{
    {21, 2705, 64, 3634401, 3635520, 3105761,   795074815,    1174, 22},
    {10, 1549, 128, 991360, 991360, 833218145, 213303845119, 1206, 20},
    {6,  960,  128, 368512, 368640, 1041300784,266573000703, 1198, 18},
    {6,  530,  128, 203520, 203520, 563129990, 144161277439, 1193, 17},
    {5,  464,  0,   148480, 148480, 442476261, 113273922816, 1184, 17}
}};
inline constexpr size_t attempt_limit = 160;
using WordSource = std::function<U()>;
using ChallengeSource = std::function<RV(size_t, bool)>;

// The verifier retains its own coins for the one submitted path at each layer.
struct LayerCoins {
  std::vector<Bytes> skipped_seeds;
  std::vector<RV> skipped_challenges;
  Bytes projection_seed;
  E zeta{};
  RV challenge;
};
struct LayerMessages {
  RV u1, u2;
  SV image;
  TerminalResponse pivot;
};

using LayerPath = std::array<Config, 5>;

// Bill verifier challenge generation separately even inside the interactive prover.
template<class F>
inline auto account(double *seconds, F operation) {
  if (!seconds) return operation();
  const auto begin = std::chrono::steady_clock::now();
  auto result = operation();
  *seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
  return result;
}

inline Bytes projection_seed(const Config &cfg, const WordSource &word) {
  const size_t bits = cfg.b + cfg.h * (3 * cfg.b - 1);
  Bytes seed((bits + 7) / 8);
  for (auto &byte : seed)
    byte = static_cast<uint8_t>(word());
  if (bits % 8)
    seed.back() &= (1u << (bits % 8)) - 1;
  return seed;
}

inline E aggregation_coin(const WordSource &word) {
  E value;
  for (auto &limb : value) {
    do { limb = word() & mask; } while (limb >= q);
  }
  return value;
}

inline Aggregated prove_projection(const Config &cfg, const System &system,
                                   const RV &w, const WordSource &word,
                                   LayerCoins &coins, LayerMessages &msg, double *coin_seconds) {
  for (size_t i = 0; i < attempt_limit; ++i) {
    coins.projection_seed = account(coin_seconds, [&]{ return projection_seed(cfg, word); });
    auto projection = expand_projection(cfg, coins.projection_seed);
    auto image = project(projection, w);
    if (proj_norm(image, cfg.S)) {
      std::cerr << "projection_attempts=" << i + 1 << '\n';
      msg.image = std::move(image);
      coins.zeta = account(coin_seconds, [&]{ return aggregation_coin(word); });
      return aggregate(system, projection, msg.image, coins.zeta, one);
    }
    coins.skipped_seeds.push_back(std::move(coins.projection_seed));
  }
  throw std::runtime_error("projection retry budget exhausted");
}

inline RV prove_response(const Config &cfg, const RV &w, const ChallengeSource &draw,
                         LayerCoins &coins, double *coin_seconds) {
  for (size_t i = 0; i < attempt_limit; ++i) {
    coins.challenge = account(coin_seconds, [&]{ return draw(cfg.s, cfg.rho == 0); });
    check(coins.challenge.size() == cfg.s, "short challenge block count");
    if (!cfg.rho)
      check(!rzero(coins.challenge.front()), "pivot challenge must be nonzero");
    auto z = response(w, coins.challenge, cfg.n);
    if (energy(z) <= cfg.G) {
      std::cerr << "short_challenge_attempts=" << i + 1 << '\n';
      return z;
    }
    coins.skipped_challenges.push_back(std::move(coins.challenge));
  }
  throw std::runtime_error("short-response retry budget exhausted");
}

struct LayerKeys {
  Key A, B, D;
  int ell;
  LayerKeys(const Config &cfg, const std::shared_ptr<const Master> &crs)
      : A{crs, 16, cfg.n}, B{crs, 6, 0}, D{crs, 6, 0},
        ell(cfg.rho ? ellrho(cfg.rho) : 0) {
    B.cols = ell * 16 * cfg.s;
    D.cols = ell * 4 * (cfg.s * (cfg.s + 1) / 2);
  }
};

inline RV child_opening(const RV &z, I radix, const RV &bt, const RV &bh) {
  RV child(2 * z.size());
  for (size_t j = 0; j < z.size(); ++j) {
    for (size_t t = 0; t < degree; ++t) {
      const I x = center(z[j][t]);
      const I magnitude = (std::llabs(x) + radix / 2 - 1) / radix;
      const I hi = x < 0 ? -magnitude : magnitude;
      child[j][t] = qm(x - radix * hi);
      child[z.size() + j][t] = qm(hi);
    }
  }
  child.insert(child.end(), bt.begin(), bt.end());
  child.insert(child.end(), bh.begin(), bh.end());
  return child;
}

inline System next_system(size_t index, const LayerKeys &keys, Aggregated ag,
                          const LayerMessages &msg, const LayerCoins &coins,
                          const LayerPath &path = layers) {
  auto child = childsystem({path[index], keys.A, keys.B, keys.D, msg.u1, msg.u2,
                            coins.challenge, std::move(ag), keys.B.cols,
                            keys.D.cols, keys.ell});
  const auto &next = path.at(index + 1);
  check(child.raw == next.raw, "profile child length");
  if (next.raw == next.N)
    return child;
  return align(child, next);
}

inline void prove_layers(System sys, RV w, const std::shared_ptr<const Master> &crs,
                         const WordSource &word, const ChallengeSource &draw,
                         std::array<LayerMessages, 5> &proof,
                         std::array<LayerCoins, 5> &tape,
                         const LayerPath &path = layers, double *coin_seconds = nullptr) {
  Key pivot_key{crs, 2, 256};
  for (size_t i = 0; i < path.size(); ++i) {
    const auto &cfg = path[i];
    auto &msg = proof[i];
    auto &coins = tape[i];
    std::cerr << "Proving backend layer " << i + 1 << '/' << path.size() << '\n';
    check(w.size() == cfg.N / degree && sys.raw == cfg.N, "aligned prover dimensions");
    check(energy(w) <= cfg.S && sys.apply(w) == sys.y, "parent witness relation");
    LayerKeys keys(cfg, crs);
    RV tags, bt;
    for (size_t j = 0; j < cfg.s; ++j) {
      auto image = keys.A.apply(slice(w, j * cfg.n, cfg.n));
      tags.insert(tags.end(), image.begin(), image.end());
    }
    if (cfg.rho) {
      bt = digits(tags, cfg.rho, keys.ell);
      msg.u1 = keys.B.apply(bt);
    } else {
      msg.pivot.pivot = pivot_key.apply(digits(slice(tags, 0, keys.A.rows), 8, 16));
      msg.pivot.tothers.assign(tags.begin() + keys.A.rows, tags.end());
    }
    auto ag = prove_projection(cfg, sys, w, word, coins, msg, coin_seconds);
    auto h = symmetric(ag, w, cfg.s, cfg.n);
    if (!cfg.rho) {
      msg.pivot.h = std::move(h);
      const size_t triangle = cfg.s * (cfg.s + 1) / 2;
      for (size_t e = 0; e < 4; ++e)
        msg.pivot.h[e * triangle][0] = 0;
      msg.pivot.z = prove_response(cfg, w, draw, coins, coin_seconds);
      return;
    }
    auto bh = digits(h, cfg.rho, keys.ell);
    msg.u2 = keys.D.apply(bh);
    auto z = prove_response(cfg, w, draw, coins, coin_seconds);
    w = child_opening(z, cfg.rho, bt, bh);
    sys = next_system(i, keys, std::move(ag), msg, coins, path);
    check(w.size() * degree == path[i + 1].raw, "child opening dimensions");
    w.resize(path[i + 1].N / degree);
  }
}

inline void check_canonical_rings(const RV &rings) {
  for (const auto &r : rings)
    check(std::all_of(r.begin(), r.end(), [](U v) { return v < q; }),
          "noncanonical ring coefficient");
}

// Replay uses the verifier's retained challenges, the public statement and the
// prover's messages. No opening or prover-computed aggregate crosses this API.
inline void verify_layers(System sys, const std::shared_ptr<const Master> &crs,
                          const std::array<LayerMessages, 5> &proof,
                          const std::array<LayerCoins, 5> &tape,
                          const LayerPath &path = layers, double *tamper_seconds = nullptr) {
  for (size_t i = 0; i < path.size(); ++i) {
    std::cerr << "Verifying backend layer " << i + 1 << '/' << path.size() << '\n';
    const auto &cfg = path[i];
    const auto &msg = proof[i];
    const auto &coins = tape[i];
    check(std::all_of(msg.image.begin(), msg.image.end(), [](U v) { return v < q; }),
          "noncanonical projection coefficient");
    check(proj_norm(msg.image, cfg.S), "projection image norm or length");
    auto ag = aggregate(sys, expand_projection(cfg, coins.projection_seed),
                        msg.image, coins.zeta, one);
    LayerKeys keys(cfg, crs);
    if (cfg.rho) {
      check(msg.pivot.pivot.empty() && msg.pivot.tothers.empty() &&
                msg.pivot.h.empty() && msg.pivot.z.empty(), "unexpected pivot fields");
      check_canonical_rings(msg.u1);
      check_canonical_rings(msg.u2);
      sys = next_system(i, keys, std::move(ag), msg, coins, path);
    } else {
      check(msg.u1.empty() && msg.u2.empty(), "unexpected child commitments");
      for (const auto *field : {&msg.pivot.pivot, &msg.pivot.tothers,
                                &msg.pivot.h, &msg.pivot.z})
        check_canonical_rings(*field);
      const Key pivot_key{crs, 2, 256};
      verifyterminal(cfg, keys.A, pivot_key, ag, msg.pivot, coins.challenge);
      const auto tamper_start = std::chrono::steady_clock::now();
      auto altered = msg.pivot;
      altered.pivot[0][0] = add(altered.pivot[0][0], 1);
      try {
        verifyterminal(cfg, keys.A, pivot_key, ag, altered, coins.challenge);
      } catch (const std::runtime_error &e) {
        if (std::string_view(e.what()) != "pivot binding equation")
          throw;
        if (tamper_seconds)
          *tamper_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-tamper_start).count();
        return;
      }
      throw std::runtime_error("verifier accepted the altered pivot commitment");
    }
  }
  throw std::runtime_error("missing terminal layer");
}

inline void run_five_layers(System aligned_root, RV witness,
                            std::shared_ptr<const Master> crs, WordSource word,
                            ChallengeSource draw_challenges) {
  check(crs && crs->data.size() == 44352, "five-layer CRS length");
  check(aligned_root.raw == layers.front().N && aligned_root.y.size() == 264036,
        "five-layer root shape");
  std::array<LayerMessages, 5> proof;
  std::array<LayerCoins, 5> tape;
  const auto start = std::chrono::steady_clock::now();
  prove_layers(aligned_root, std::move(witness), crs, word, draw_challenges, proof, tape);
  const auto proved = std::chrono::steady_clock::now();
  verify_layers(std::move(aligned_root), crs, proof, tape);
  std::cerr << "backend_prove_seconds=" << std::chrono::duration<double>(proved-start).count()
            << " backend_verify_and_tamper_seconds="
            << std::chrono::duration<double>(std::chrono::steady_clock::now()-proved).count() << '\n';
}

} // namespace sonata::example
