// Independent local checks of the child rows and terminal pivot opening.
#include "fixture.hpp"
#include "sonata/backend.hpp"
#include <iostream>

using namespace sonata;
using namespace sonata::fixture;

int main() {
  constexpr size_t s = 2, n = 1, rho = 64;
  const int ell = ellrho(rho);
  const size_t bt = ell * 16 * s, bh = ell * 4 * (s * (s + 1) / 2);
  const Config cfg{s, n, rho, 64 * s * n, 64 * s * n, 100, 10000, 0, 0};
  auto master = std::make_shared<Master>();
  std::mt19937_64 rng(92);
  master->data.resize(6 * std::max(bt, bh));
  for (auto &cell : master->data)
    cell = randR(rng);
  Key A{master, 16, n}, B{master, 6, bt}, D{master, 6, bh},
      pivot{master, 2, 256};
  RV w(s * n), c(s);
  w[0][0] = 2;
  w[1][1] = 3;
  c[0][0] = 1;
  c[1][1] = 1;
  Aggregated ag;
  for (size_t e = 0; e < 4; ++e) {
    ag.phi[e].resize(s * n);
    for (auto &beta : ag.phi[e])
      beta = randR(rng);
    ag.rhs[e] = ctDot(ag.phi[e], w);
  }
  const RV h = symmetric(ag, w, s, n);
  const RV z = response(w, c, n);
  const RV t0 = A.apply({w[0]}), t1 = A.apply({w[1]});
  RV t = t0;
  t.insert(t.end(), t1.begin(), t1.end());
  const auto t_digits = digits(t, rho, ell), h_digits = digits(h, rho, ell);
  const auto child = childsystem(
      {cfg, A, B, D, B.apply(t_digits), D.apply(h_digits), c, ag, bt, bh, ell});
  RV opening = z;
  opening.resize(2 * n); // high part is zero for this short response
  opening.insert(opening.end(), t_digits.begin(), t_digits.end());
  opening.insert(opening.end(), h_digits.begin(), h_digits.end());
  check(child.apply(opening) == child.y, "child witness");
  opening[2 * n][0] = add(opening[2 * n][0], 1);
  check(child.apply(opening) != child.y, "modified child opening");
  checkadj(child, opening, rng);

  TerminalResponse last{t1, pivot.apply(digits(t0, 8, 16)), h, z};
  for (size_t e = 0; e < 4; ++e)
    last.h[e * s * (s + 1) / 2][0] = 0;
  auto terminal_cfg = cfg;
  terminal_cfg.rho = 0;
  verifyterminal(terminal_cfg, A, pivot, ag, last, c);
  last.pivot[0][0] = add(last.pivot[0][0], 1);
  bool rejected = false;
  try {
    verifyterminal(terminal_cfg, A, pivot, ag, last, c);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  check(rejected, "modified pivot commitment");
  std::cout
      << "PASS child rows and adjoint, pivot acceptance and mutation checks\n";
}
