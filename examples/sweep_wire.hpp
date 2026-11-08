#pragma once
#include "sonata/backend.hpp"
#include <limits>

namespace sonata::example {
namespace wire_detail {

// Lexicographic (most significant limb first). The ambient ball has <2^160 points.
using Rank = std::array<U, 3>;
inline void add_rank(Rank &a, const Rank &b) {
  U carry = 0;
  for (int i = 2; i >= 0; --i) {
    const W v = W(a[i]) + b[i] + carry;
    a[i] = U(v);
    carry = U(v >> 64);
  }
  check(!carry, "ball count overflow");
}
inline void subtract_rank(Rank &a, const Rank &b) {
  U borrow = 0;
  for (int i = 2; i >= 0; --i) {
    const W subtrahend = W(b[i]) + borrow;
    borrow = W(a[i]) < subtrahend;
    a[i] -= U(subtrahend);
  }
  check(!borrow, "ball rank underflow");
}
inline int radius(int budget) {
  int v = 0;
  while ((v + 1) * (v + 1) <= budget) ++v;
  return v;
}
struct Ball {
  std::array<std::array<Rank, 129>, 65> counts{};
  Ball() {
    counts[0].fill(Rank{0, 0, 1});
    for (size_t n = 1; n <= 64; ++n)
      for (int b = 0; b <= 128; ++b) {
        const int r = radius(b);
        for (int x = -r; x <= r; ++x)
          add_rank(counts[n][b], counts[n-1][b-x*x]);
      }
  }
  Rank rank(const Ring &c) const {
    Rank result{};
    int b = 128;
    for (int j = 0; j < 64; ++j) {
      check(c[j] < q, "noncanonical challenge");
      const I v = center(c[j]);
      const int r = radius(b);
      check(v >= -r && v <= r, "challenge outside ball");
      for (int x = -r; x < v; ++x) add_rank(result, counts[63-j][b-x*x]);
      b -= int(v*v);
    }
    return result;
  }
  Ring unrank(Rank rank) const {
    check(rank < counts[64][128], "rank outside ambient ball");
    Ring c{};
    int b = 128;
    for (int j = 0; j < 64; ++j) {
      const int r = radius(b);
      for (int x = -r; x <= r; ++x) {
        const auto &weight = counts[63-j][b-x*x];
        if (rank < weight) {
          c[j] = qm(x);
          b -= x*x;
          break;
        }
        subtract_rank(rank, weight);
      }
    }
    return c;
  }
};
inline const Ball& ball() {
  static const Ball index;
  return index;
}

struct Rice {
  size_t bytes;
  unsigned h;
};
inline Rice rice_layout(size_t n, U energy_limit) {
  check(n > 0 && n <= 10000000 && energy_limit < (U(1)<<62), "bounded vector shape");
  U best_bits = std::numeric_limits<U>::max();
  unsigned best_h = 1;
  U budget = energy_limit;
  for (unsigned h = 1;; ++h, budget /= 4) {
    // The cheapest way for n nonnegative integer quotients to sum to t
    // distributes t as evenly as possible. Zigzag(x)>>h <= |x|/2^(h-1).
    U lo = 0, hi = budget+1;
    while (hi-lo > 1) {
      const U t = lo + (hi-lo)/2, a = t/n, rem = t%n;
      const W minimum = W(t)*a + W(rem)*(a+1);
      if (minimum <= budget) lo=t; else hi=t;
    }
    const U bits = n*(h+1)+lo;
    if (bits < best_bits) { best_bits=bits; best_h=h; }
    if (!budget) break;
  }
  return {size_t((best_bits+7)/8), best_h};
}
} // namespace wire_detail

class WireWriter {
  void scalar(U v) {
    check(v<q, "noncanonical field coefficient");
    for (unsigned j=0; j<6; ++j) data.push_back(uint8_t(v>>(8*j)));
  }
public:
  Bytes data;
  void field(const E &v) { for (U a:v) scalar(a); }
  void rings(const RV &v) { for (const auto &r:v) for (U a:r) scalar(a); }
  void seed(const Bytes &v) { data.insert(data.end(),v.begin(),v.end()); }
  void flag(bool submit) { data.push_back(submit ? 1 : 0); }
  void triangle(const RV &v, size_t s) {
    const size_t tri=s*(s+1)/2;
    check(s>0 && v.size()==4*tri, "terminal triangle shape");
    for (size_t j=0; j<v.size(); ++j)
      for (size_t i=0; i<degree; ++i) {
        if (j%tri || i) scalar(v[j][i]);
        else check(v[j][i]==0, "terminal placeholder must be zero");
      }
  }
  void challenges(const RV &v) {
    for (const auto &r:v) {
      auto rank=wire_detail::ball().rank(r);
      for (size_t i=0; i<20; ++i) data.push_back(uint8_t(rank[2-i/8]>>(8*(i%8))));
    }
  }
  void bounded(const SV &v, U energy_limit) {
    const auto layout=wire_detail::rice_layout(v.size(),energy_limit);
    size_t pos=8*data.size();
    data.resize(data.size()+layout.bytes,0);
    const size_t end=8*data.size();
    W used=0;
    for (U a:v) {
      check(a<q,"noncanonical bounded coefficient");
      const I x=center(a);
      const U mag=U(x<0 ? -x : x);
      used+=W(mag)*mag;
      check(used<=energy_limit,"bounded vector energy");
      const U z=2*mag-U(x<0);
      const U zeros=z>>layout.h;
      check(zeros < end-pos && layout.h<end-pos-zeros,"Rice output exceeds padding");
      pos+=size_t(zeros);
      data[pos/8]|=uint8_t(1u<<(pos%8));
      ++pos;
      for (unsigned j=0; j<layout.h; ++j,++pos)
        data[pos/8]|=uint8_t(((z>>j)&1)<<(pos%8));
    }
  }
};

class WireReader {
  const Bytes &bytes_;
  size_t pos_=0;
  void need(size_t n) const { check(n<=bytes_.size()-pos_,"truncated wire message"); }
  U scalar() {
    need(6);
    U v=0;
    for (unsigned j=0; j<6; ++j) v|=U(bytes_[pos_++])<<(8*j);
    check(v<q,"noncanonical wire field");
    return v;
  }
public:
  explicit WireReader(const Bytes &bytes):bytes_(bytes) {}
  bool done() const { return pos_==bytes_.size(); }
  E field() { E v; for (auto &a:v) a=scalar(); return v; }
  RV rings(size_t n) {
    check(n<=(bytes_.size()-pos_)/(degree*6),"truncated wire rings");
    RV v(n); for (auto &r:v) for (auto &a:r) a=scalar(); return v;
  }
  Bytes seed(size_t n) {
    need(n); Bytes out(bytes_.begin()+pos_,bytes_.begin()+pos_+n); pos_+=n; return out;
  }
  bool flag() { need(1); auto x=bytes_[pos_++]; check(x<=1,"unknown stage control"); return x==1; }
  RV triangle(size_t s) {
    const size_t tri=s*(s+1)/2;
    check(s>0,"invalid terminal arity");
    need((4*tri*degree-4)*6);
    RV v(4*tri);
    for (size_t j=0; j<v.size(); ++j)
      for (size_t i=0; i<degree; ++i)
        if (j%tri || i) v[j][i]=scalar();
    return v;
  }
  RV challenges(size_t n) {
    check(n<=(bytes_.size()-pos_)/20,"truncated wire challenges");
    RV v(n);
    for (auto &r:v) {
      wire_detail::Rank rank{};
      for (size_t i=0; i<20; ++i) rank[2-i/8]|=U(bytes_[pos_++])<<(8*(i%8));
      r=wire_detail::ball().unrank(rank);
    }
    return v;
  }
  SV bounded(size_t n, U energy_limit) {
    const auto layout=wire_detail::rice_layout(n,energy_limit);
    need(layout.bytes);
    size_t bit=8*pos_, end=bit+8*layout.bytes;
    auto read_bit=[&]() -> U {
      check(bit<end,"truncated Rice payload");
      const size_t at=bit++;
      return (bytes_[at/8]>>(at%8))&1;
    };
    SV v(n);
    W used=0;
    for (U &a:v) {
      U zeros=0;
      while (!read_bit()) ++zeros;
      check(zeros<=q>>layout.h,"oversized Rice quotient");
      U z=zeros<<layout.h;
      for (unsigned j=0;j<layout.h;++j) z|=read_bit()<<j;
      const U mag=(z+1)/2;
      check(mag<=q/2,"noncanonical signed Rice value");
      used+=W(mag)*mag;
      check(used<=energy_limit,"decoded vector exceeds energy");
      a=qm((z&1) ? -I(mag) : I(mag));
    }
    while (bit<end) check(!read_bit(),"nonzero Rice padding");
    pos_+=layout.bytes;
    return v;
  }
};
} // namespace sonata::example
