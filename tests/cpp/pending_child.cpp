// Algebra regression for the pending-root to nonterminal-child boundary.
#include "fixture.hpp"
#include "sonata/backend.hpp"
#include <iostream>

using namespace sonata;
using namespace sonata::fixture;

static PublicFront committed_transition(const Fixture &f,
                                        const std::shared_ptr<Master> &master) {
  PublicFront pub{};
  pub.k = f.k;
  pub.alpha = f.weights.front();
  pub.rho = f.rho;
  pub.gamma = f.gamma;
  pub.eta = f.eta;
  pub.prev_pending = {9, 0, 0, 0};
  pub.pending = ea(pub.prev_pending, f.pending);
  check(pub.pending != pub.prev_pending, "expected a nonzero pending update");
  Key state{master, 4, 219};
  for (const auto &s : f.states) {
    std::array<E, 5> fields;
    std::copy_n(s.begin(), fields.size(), fields.begin());
    pub.pub.push_back(fields);
    RV opening(state.cols);
    for (size_t a = 0; a < s.size(); ++a)
      for (size_t t = 0; t < 4; ++t)
        for (size_t j = 0; j < 48; ++j)
          put(opening, 192*a + 48*t + j, (s[a][t] >> j) & 1);
    pub.cs.push_back(state.apply(opening));
  }
  const size_t groups = f.scalarvalues.size();
  Key aux{master, 4, (7*groups + 63)/64};
  RV auxiliary(aux.cols);
  for (size_t g = 0; g < groups; ++g)
    for (size_t j = 0; j < 7; ++j)
      put(auxiliary, 7*g + j, f.bits[55*g + 48 + j]);
  pub.caux = aux.apply(auxiliary);
  E claim = ea(E{mul(groups, KAPPA), 0, 0, 0}, em(pub.eta, f.pending));
  pub.field = prove(f.table(), claim, 75011, true);
  return pub;
}

// Deterministic verifier seed material for the archived k=4 dimensions.
static Projection fixture_projection(const Config &cfg, uint64_t fixture_seed) {
  std::mt19937_64 rng(fixture_seed);
  const size_t bits = cfg.b + cfg.h * (3 * cfg.b - 1);
  Bytes seed((bits + 7) / 8);
  for (auto &v : seed)
    v = static_cast<uint8_t>(rng());
  if (bits % 8)
    seed.back() &= (1u << (bits % 8)) - 1;
  return expand_projection(cfg, seed);
}

static void check_pivot(const Config &cfg, const System &child, RV witness,
                        const std::shared_ptr<Master> &master) {
  const auto system = align(child, cfg);
  witness.resize(cfg.N/64); // Padding rows authenticate zeros.
  check(system.apply(witness) == system.y && energy(witness) <= cfg.S,
        "same nonterminal child opening at pivot input");
  Key A{master, 16, cfg.n}, pivot{master, 2, 256};
  const auto t0 = A.apply(slice(witness, 0, cfg.n));
  TerminalResponse terminal;
  terminal.pivot = pivot.apply(digits(t0, 8, 16));
  for (size_t j = 1; j < cfg.s; ++j) {
    auto block = A.apply(slice(witness, j*cfg.n, cfg.n));
    terminal.tothers.insert(terminal.tothers.end(), block.begin(), block.end());
  }
  const auto projection = fixture_projection(cfg, 8702);
  const auto image = project(projection, witness);
  check(proj_norm(image, cfg.S), "pivot projection bound");
  const auto ag = aggregate(system, projection, image, E{7,2,3,1}, one);
  terminal.h = symmetric(ag, witness, cfg.s, cfg.n);
  const size_t ht = cfg.s*(cfg.s+1)/2;
  for (size_t e = 0; e < 4; ++e)
    terminal.h[e*ht][0] = 0;
  RV c(cfg.s);
  c[0][0] = 1;
  c[1][1] = 1;
  c[2][0] = q-1;
  terminal.z = response(witness, c, cfg.n);
  verifyterminal(cfg, A, pivot, ag, terminal, c);
  terminal.pivot[0][0] = add(terminal.pivot[0][0], 1);
  bool rejected = false;
  try {
    verifyterminal(cfg, A, pivot, ag, terminal, c);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected, "altered public pivot commitment");
}

int main() {
  // The archived k=4 first layer has three blocks and radix 64.
  const Config cfg{3, 230, 64, 44001, 44160, 37601, 9625855, 1081, 15};
  const Config terminal_cfg{3, 346, 0, 66304, 66432, 52543624,
                            13451167743, 1149, 16};
  const size_t n = cfg.n, s = cfg.s, ht = s * (s+1)/2;
  const int ell = ellrho(cfg.rho);
  std::mt19937_64 rng(707070);
  auto master = std::make_shared<Master>();
  const size_t bt = ell*16*s, bh = ell*4*ht;
  master->data.resize(std::max({16*terminal_cfg.n, 16*n, 6*bt, 6*bh,
                                size_t(4*219), size_t(4*((7*(64*4+544)+63)/64))}));
  for (auto &cell : master->data)
    cell = randR(rng);
  Fixture instance(4, 511, true);
  auto pub = committed_transition(instance, master);
  const auto raw = initialsystem(pub, master);
  const auto root = align(raw, cfg);
  RV w(cfg.N/64);
  for (size_t i = 0; i < instance.bits.size(); ++i)
    put(w, i, instance.bits[i]);
  put(w, raw.raw-1, 1);
  check(root.apply(w) == root.y, "original pending handoff");

  Key A{master, 16, n}, B{master, 6, bt}, D{master, 6, bh};
  RV t;
  for (size_t j = 0; j < s; ++j) {
    auto block_image = A.apply(slice(w, j*n, n));
    t.insert(t.end(), block_image.begin(), block_image.end());
  }
  auto t_digits = digits(t, cfg.rho, ell);
  const auto u1 = B.apply(t_digits);
  const auto projection = fixture_projection(cfg, 8701);
  auto image = project(projection, w);
  check(proj_norm(image, cfg.S), "projection image");
  const auto ag = aggregate(root, projection, image, E{1,5,2,3}, one);
  auto h = symmetric(ag, w, s, n);
  const auto h_original = digits(h, cfg.rho, ell);
  const auto u2 = D.apply(h_original);
  RV c(s);
  c[0][0] = 1;
  c[1][1] = 1;
  c[2][0] = q-1;
  auto z = response(w, c, n);
  for (const auto &ring : z)
    for (U coeff : ring)
      check(std::abs(center(coeff)) <= 3, "fixture response fits low digit");
  RV prefix = z;
  prefix.resize(2*n); // Every response fits in its low digit.
  prefix.insert(prefix.end(), t_digits.begin(), t_digits.end());
  const size_t principal_offset = 64*(B.rows + D.rows + A.rows);
  check(principal_offset == 1792, "public row layout");
  for (size_t lane = 0; lane <= 4; ++lane) {
    auto table = h;
    if (lane < 4) {
      auto &entry = table[lane*ht + tri(s,0,1)][0];
      entry = add(entry, 1);
    }
    // Recommit the perturbed off-diagonal so commitment checks still pass.
    const auto h_digits = lane == 4 ? h_original : digits(table, cfg.rho, ell);
    const auto tag = lane == 4 ? u2 : D.apply(h_digits);
    const auto child = childsystem({cfg, A, B, D, u1, tag,
                                    c, ag, bt, bh, ell});
    auto opening = prefix;
    opening.insert(opening.end(), h_digits.begin(), h_digits.end());
    auto expected = child.y;
    if (lane < 4) {
      // Principal LHS minus RHS changes by -2*c0*c1 = -2X only here.
      const size_t row = principal_offset + 64*lane + 1;
      expected[row] = sub(expected[row], 2);
    }
    check(child.apply(opening) == expected, "pending child row residuals");
    if (lane == 4) {
      check_pivot(terminal_cfg, child, opening, master);
      for (auto &ring : opening)
        ring = randR(rng);
      checkadj(child, opening, rng);
    }
  }
  std::cout << "PASS pending root to four-lane child and its pivot, exact row mutations and adjoint\n";
}
