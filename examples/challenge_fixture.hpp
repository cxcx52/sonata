#pragma once

#include "sonata/ring.hpp"
#include <iostream>

namespace sonata::example {

// The Python verifier samples from the spectral support. The pipe carries
// signed coefficients of its challenge. The small live run uses OS randomness.
inline RV draw_challenges(size_t blocks, bool terminal) {
  std::cout << "SONATA_DRAW " << blocks << ' ' << int(terminal) << std::endl;
  RV challenges(blocks);
  for (auto &ring : challenges) {
    unsigned squared_norm = 0;
    for (auto &coefficient : ring) {
      int value;
      check(bool(std::cin >> value), "incomplete challenge from verifier");
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
