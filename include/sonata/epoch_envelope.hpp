#pragma once
#include "sonata/backend.hpp"

namespace sonata::epoch {

// The active relation takes ceil(active.raw/64) rings and may share its last
// ring with ignored coordinates. The imported align() still enforces its own
// final N - raw padding rows. Pass the result of this function to align().
inline System embed_prefix(System active, const Config &cfg, size_t aligned_rows) {
  check(active.raw > 0 && active.raw <= cfg.raw && cfg.raw <= cfg.N &&
            cfg.N % 64 == 0, "epoch coefficient capacity");
  const size_t fixed_padding = cfg.N - cfg.raw;
  check(aligned_rows >= fixed_padding &&
            active.y.size() <= aligned_rows - fixed_padding, "epoch row capacity");
  const size_t rows = aligned_rows - fixed_padding;
  const size_t rings = (active.raw + 63) / 64;
  const size_t last = active.raw % 64;
  const auto source = std::make_shared<const System>(std::move(active));
  System prefix;
  prefix.raw = cfg.raw;
  prefix.y = source->y;
  prefix.y.resize(rows, 0);
  prefix.apply = [source, rows, rings, last](const RV &w) {
    check(w.size() >= rings, "short epoch input");
    RV v(w.begin(), w.begin() + rings);
    if (last) std::fill(v.back().begin() + last, v.back().end(), U(0));
    auto out = source->apply(v);
    check(out.size() == source->y.size(), "active epoch output dimension");
    out.resize(rows, 0);
    return out;
  };
  prefix.adj = [source, rows, rings, last, raw = cfg.raw](const SV &a) {
    check(a.size() == rows, "epoch adjoint row capacity");
    auto out = source->adj(SV(a.begin(), a.begin() + source->y.size()));
    check(out.size() == rings, "active epoch adjoint dimension");
    // Ring rows are already adjointed for the constant-term pairing.
    if (last) {
      auto coefficients = dual(out.back());
      std::fill(coefficients.begin() + last, coefficients.end(), U(0));
      out.back() = dual(coefficients);
    }
    out.resize((raw + 63) / 64);
    return out;
  };
  return prefix;
}

} // namespace sonata::epoch
