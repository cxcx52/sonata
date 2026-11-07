#pragma once
#include "sonata/field.hpp"
#include <limits>

namespace sonata::epoch {

// In adjacent-cell folding, r[0] is the least significant address bit.
// Equality weights factor across a split in this address. The product uses
// no inverses, so the same formula applies to verifier coordinates 0 and 1.
class FactoredWeights {
  std::vector<E> low_, high_;
  size_t low_bits_;
public:
  explicit FactoredWeights(const std::vector<E>& r) : low_bits_(r.size()/2) {
    check(r.size()<std::numeric_limits<size_t>::digits,"weight dimension");
    const auto mid=r.begin()+low_bits_;
    low_=eqweights(std::vector<E>(r.begin(),mid));
    high_=eqweights(std::vector<E>(mid,r.end()));
  }
  size_t size() const { return low_.size()*high_.size(); }
  E operator[](size_t address) const {
    return em(low_[address & (low_.size()-1)], high_.at(address >> low_bits_));
  }
};
} // namespace sonata::epoch
