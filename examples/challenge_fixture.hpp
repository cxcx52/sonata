#pragma once

#include "sonata/ring.hpp"
#include <iostream>

namespace sonata::example {

// run_e2e.py owns the deterministic verifier coins and checks spectral membership
// with the reference sampler. The pipe carries signed coefficients for this demo.
inline RV draw_challenges(size_t blocks, bool terminal) {
  std::cout << "SONATA_DRAW " << blocks << ' ' << int(terminal) << std::endl;
  RV challenges(blocks);
  for (auto &ring : challenges) {
    unsigned squared_norm = 0;
    for (auto &coefficient : ring) {
      int value;
      check(bool(std::cin >> value), "incomplete challenge from run_e2e.py");
      check(value >= -11 && value <= 11, "short-challenge coefficient");
      squared_norm += value * value;
      coefficient = qm(value);
    }
    check(squared_norm <= 128, "short-challenge energy");
  }
  if (terminal)
    check(!challenges.empty() && !rzero(challenges.front()), "zero pivot challenge");
  return challenges;
}

} // namespace sonata::example
