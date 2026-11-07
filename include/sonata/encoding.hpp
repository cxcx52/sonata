#pragma once

#include "field.hpp"

namespace sonata {

// The first 48 coordinates encode the field representative. In the auxiliary
// mask, ties on Hamming weight take the smaller integer mask.
inline std::array<uint8_t, 55> canonical_group(U value) {
  check(value < q, "noncanonical input representative");
  std::array<uint8_t, 55> bits;
  for (size_t j = 0; j < 48; ++j)
    bits[j] = (value >> j) & 1;
  const int zeros = 42 - __builtin_popcountll(value >> 6);
  const int exceptional = zeros == 0;
  const int target = (zeros - 1 + exceptional) / 2;
  int best = -1;
  for (int candidate = 0; candidate < 32; ++candidate) {
    int weight = 0;
    for (int j = 0; j < 5; ++j)
      if ((candidate >> j) & 1)
        weight += HWEIGHTS[j];
    if (weight == target &&
        (best < 0 || __builtin_popcount(unsigned(candidate)) <
                         __builtin_popcount(unsigned(best))))
      best = candidate;
  }
  check(best >= 0, "aux count representability");
  for (int j = 0; j < 5; ++j)
    bits[48 + j] = (best >> j) & 1;
  bits[53] = exceptional;
  bits[54] = exceptional * bits[2];
  return bits;
}

} // namespace sonata
