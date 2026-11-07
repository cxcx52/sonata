#include "fixture.hpp"
#include "sonata/backend.hpp"

#include <chrono>
#include <iostream>
#include <string>

using namespace sonata;
using Clock = std::chrono::steady_clock;

static auto elapsed(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
}

static void encoding() {
  constexpr size_t groups = 66064;
  std::mt19937_64 rng(7931);
  SV inputs(groups);
  for (auto& x : inputs) x = fixture::randomQ(rng);
  inputs.front() = 0;
  inputs.back() = q - 1;
  RV opening((55 * groups + 63) / 64);
  const auto start = Clock::now();
  for (size_t i = 0; i < groups; ++i) {
    const auto encoded = canonical_group(inputs[i]);
    for (size_t j = 0; j < encoded.size(); ++j)
      put(opening, i * 55 + j, encoded[j]);
  }
  const auto stop = Clock::now();
  for (size_t g = 0; g < groups; ++g) {
    U decoded = 0;
    for (size_t j = 0; j < 48; ++j)
      decoded |= cv(opening, 55 * g + j) << j;
    check(decoded == inputs[g], "encoding round trip");
  }
  std::cout << "{\"wall_ns\":" << elapsed(start, stop) << ",\"scalars\":" << groups
            << ",\"coefficients\":" << 55 * groups << ",\"energy\":" << U(energy(opening)) << "}\n";
}

static void sumcheck() {
  // Materialize an honest four-source square table before starting the clock.
  fixture::Fixture f(4, 511);
  auto table = f.table();
  const size_t entries = table.size();
  const E claim{mul(f.scalarvalues.size(), KAPPA), 0, 0, 0};
  const auto start = Clock::now();
  const auto proof = fixture::prove(std::move(table), claim, 75011, true);
  const auto stop = Clock::now();
  check(verify(proof, claim, proof.terminal), "square Sumcheck recurrence");
  std::cout << "{\"wall_ns\":" << elapsed(start, stop) << ",\"table_entries\":" << entries
            << ",\"rounds\":" << proof.ch.size() << ",\"terminal_limb\":" << proof.terminal[0] << "}\n";
}

static void projection() {
  const Config cfg{3, 230, 64, 44001, 44160, 37601, 9625855, 1081, 15};
  std::mt19937_64 rng(7061);
  const size_t bits = cfg.b + cfg.h * (3 * cfg.b - 1);
  Bytes seed((bits + 7) / 8);
  for (auto& v : seed) v = static_cast<uint8_t>(rng());
  if (bits % 8) seed.back() &= (1u << (bits % 8)) - 1;
  RV w(cfg.N / 64);
  for (auto& ring : w)
    for (auto& x : ring) x = rng() & 1;
  const auto start = Clock::now();
  const auto p = expand_projection(cfg, seed);
  const auto middle = Clock::now();
  const auto image = project(p, w);
  const auto stop = Clock::now();
  I row0 = 0;
  for (size_t j = 0; j < cfg.N; ++j)
    row0 += p.pi[j] * center(cv(w, j));
  check(image.size() == projrows && image[0] == qm(row0), "projected first row");
  std::cout << "{\"wall_ns\":" << elapsed(start, stop)
            << ",\"expand_ns\":" << elapsed(start, middle) << ",\"apply_ns\":" << elapsed(middle, stop)
            << ",\"rows\":" << projrows << ",\"columns\":" << cfg.N << ",\"seed_bits\":" << bits
            << ",\"hash_calls\":" << p.hashcalls << ",\"row0\":" << row0 << "}\n";
}

int main(int argc, char** argv) {
  try {
    const std::string name = argc == 2 ? argv[1] : "";
    if (name == "encoding") encoding();
    else if (name == "sumcheck") sumcheck();
    else if (name == "projection") projection();
    else throw std::runtime_error("expected encoding, sumcheck or projection");
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
