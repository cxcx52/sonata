#include "dense_handoff.hpp"
#include "fixture.hpp"
#include "sonata/backend.hpp"
#include <iostream>
using namespace sonata;
using namespace sonata::fixture;

int main() {
  constexpr int arity = 2;
  const Config cfg{3, 193, 0, 36961, 37056, 31585, 8085760, 1081, 15};
  auto master = std::make_shared<Master>();
  std::mt19937_64 setup_rng(707070);
  master->data.resize(16 * cfg.n);
  for (auto &cell : master->data)
    cell = randR(setup_rng);
  Fixture instance(arity, 511, true);
  PublicFront f{};
  f.k = arity;
  f.alpha = {U(f.k + 17), 3, 4, 5};
  f.rho = instance.rho;
  f.gamma = instance.gamma;
  f.eta = instance.eta;
  f.prev_pending = {9, 0, 0, 0};
  f.pending = ea(f.prev_pending, instance.pending);
  check(f.pending != f.prev_pending,
        "fixture needs a nonzero fresh residual sketch");
  Key state{master, 4, 219};
  for (const auto &s : instance.states) {
    std::array<E, 5> pub;
    std::copy_n(s.begin(), 5, pub.begin());
    f.pub.push_back(pub);
    RV opening(219);
    for (size_t a = 0; a < 73; ++a)
      for (size_t t = 0; t < 4; ++t)
        for (size_t j = 0; j < 48; ++j)
          put(opening, 192 * a + 48 * t + j, (s[a][t] >> j) & 1);
    f.cs.push_back(state.apply(opening));
  }
  size_t V = instance.scalarvalues.size();
  Key aux{master, 4, (7 * V + 63) / 64};
  RV aux_opening(aux.cols);
  for (size_t g = 0; g < V; ++g)
    for (size_t j = 0; j < 7; ++j)
      put(aux_opening, 7 * g + j, instance.bits[55 * g + 48 + j]);
  f.caux = aux.apply(aux_opening);
  E claim =
      ea(E{mul(V, KAPPA), 0, 0, 0}, em(f.eta, es(f.pending, f.prev_pending)));
  f.field = prove(instance.table(), claim, 75011, true);
  auto raw = initialsystem(f, master);
  auto sys = align(raw, cfg);
  RV witness(cfg.N / 64);
  for (size_t j = 0; j < instance.bits.size(); ++j)
    put(witness, j, instance.bits[j]);
  put(witness, raw.raw - 1, 1);
  check(sys.apply(witness) == sys.y,
        "pending endpoint or state/auxiliary links");

  // Check the transposed public system off the honest witness as well.
  std::mt19937_64 rng(19);
  for (size_t j = 0; j < witness.size(); ++j)
    witness[j] = randR(rng);
  checkadj(sys, witness, rng);
  E endpoint = dense_handoff::check_table(f, witness, f.field.ch);
  dense_handoff::check_root(f, raw, state, aux, witness, endpoint);

  // Independent batching coins and an independently drawn multilinear point.
  // fieldadj's linear identity applies even away from a valid transcript.
  auto random_front = f;
  random_front.rho = randomE(rng);
  random_front.eta = randomE(rng);
  random_front.gamma = randomE(rng);
  random_front.alpha = randomE(rng);
  auto point = f.field.ch;
  for (auto &v : point)
    v = randomE(rng);
  dense_handoff::check_table(random_front, witness, point);
  // Zero projection isolates the L,y part of the real four-lane aggregate.
  Projection pr{};
  pr.N = cfg.N;
  pr.pi.resize(projrows * cfg.N);
  E zeta = {1, 5, 2, 3}, basis_power = one, expected = zero;
  auto ag = aggregate(sys, pr, SV(projrows), zeta, one);
  auto rows = sys.apply(witness);
  for (size_t j = 0; j < rows.size(); j += 4) {
    E packed{};
    for (size_t lane = 0; lane < 4 && j + lane < rows.size(); ++lane)
      packed[lane] = sub(rows[j + lane], sys.y[j + lane]);
    expected = ea(expected, em(basis_power, packed));
    basis_power = em(basis_power, zeta);
  }
  for (size_t e = 0; e < 4; ++e)
    check(sub(ctDot(ag.phi[e], witness), ag.rhs[e]) == expected[e],
          "off-witness four-lane residual");
  std::cout << "PASS pending endpoint, original commitments, off-witness "
               "adjoint, dense semantic table (all 36960 coefficients and "
               "offset at two points), and four lanes\n";
}
