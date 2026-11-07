#pragma once
#include "sonata/frontend.hpp"
#include <utility>

namespace dense_handoff {
using namespace sonata;

// Physical layout from the relation definition: 272 old base-field limbs,
// 64 limbs per fresh input, then 272 output limbs.
inline int group(int k, int s, int slot, int limb) {
  if (slot < 5)
    return -1;
  if (s == 0)
    return 4 * (slot - 5) + limb;
  if (s == k + 1)
    return 272 + 64 * k + 4 * (slot - 5) + limb;
  if (slot >= 69 || limb != 0)
    return -1;
  return 272 + 64 * (s - 1) + (slot - 5);
}

inline E power(E x, unsigned n) {
  E v = one;
  while (n) {
    if (n & 1U)
      v = em(v, x);
    x = sq(x);
    n >>= 1;
  }
  return v;
}

// Each table cell is a forward affine expression in the unconstrained Fq
// witness coordinates. Terms can repeat. No assumption of Boolean inputs.
struct Cell {
  E constant{};
  std::vector<std::pair<size_t, E>> terms;
  void add_slot(const PublicFront &f, int s, int slot, E factor) {
    if (slot < 5) {
      constant = ea(constant, em(factor, f.pub[s][slot]));
      return;
    }
    for (int t = 0; t < 4; ++t) {
      int g = group(f.k, s, slot, t);
      if (g < 0)
        continue;
      E basis{};
      basis[t] = 1;
      E coeff = em(basis, factor);
      for (int j = 0; j < 48; ++j)
        terms.emplace_back(55 * size_t(g) + j, sc(coeff, U(1) << j));
    }
  }
};

// Build one cell at a time. Accumulate its symbolic coefficients, and
// independently evaluate the concrete table for the multilinear folding check.
struct Oracle {
  E offset{};
  std::vector<E> row, values;
  const RV &input;
  const std::vector<E> &point;
  Oracle(size_t groups, const RV &x, const std::vector<E> &r)
      : row(55 * groups), input(x), point(r) {}
  void emit(const Cell &cell, E scale = one) {
    const size_t idx = values.size();
    E weight = one;
    for (size_t j = 0; j < point.size(); ++j)
      weight = em(weight, ((idx >> j) & 1U) ? point[j] : es(one, point[j]));
    E factor = em(weight, scale);
    offset = ea(offset, em(factor, cell.constant));
    E value = cell.constant;
    for (const auto &[col, coeff] : cell.terms) {
      row[col] = ea(row[col], em(factor, coeff));
      value = ea(value, sc(coeff, cv(input, col)));
    }
    values.push_back(em(scale, value));
  }
  E fold_table() const {
    auto table = values;
    table.resize(size_t(1) << point.size());
    for (auto r : point) {
      for (size_t j = 0; j < table.size() / 2; ++j)
        table[j] = ea(table[2 * j], em(r, es(table[2 * j + 1], table[2 * j])));
      table.resize(table.size() / 2);
    }
    return table.at(0);
  }
};

inline E check_table(const PublicFront &f, const RV &x,
                     const std::vector<E> &r) {
  const int V = 64 * f.k + 544;
  Oracle oracle(V, x, r);
  for (int g = 0; g < V; ++g)
    for (int i = 0; i < 55; ++i) {
      Cell cell;
      cell.constant[0] = BETA[i];
      // Use every actual LDL coefficient, including all high-bit columns and
      // the seven auxiliaries. The production code optimizes this with a
      // prefix.
      for (int j = 0; j < 55; ++j)
        if (LDL_L[j][i])
          cell.terms.emplace_back(55 * g + j, E{LDL_L[j][i], 0, 0, 0});
      oracle.emit(cell,
                  E{SQRT_D[i][0], SQRT_D[i][1], SQRT_D[i][2], SQRT_D[i][3]});
    }
  Cell linear;
  const U half = (q + 1) / 2;
  for (int s = 0; s < f.k + 2; ++s)
    for (int j = 0; j < 4; ++j) {
      bool fresh = s > 0 && s <= f.k;
      E w = em(f.eta, fresh ? power(f.gamma, 4 * (s - 1) + j + 1)
                            : power(f.rho, (s == 0 ? 0 : 4) + j + 1));
      for (int l = 0; l < 8; ++l) {
        int a = 5 + 16 * j + 2 * l;
        for (int sign : {1, -1}) {
          Cell cell;
          cell.add_slot(f, s, a, w);
          cell.add_slot(f, s, a + 1, sign == 1 ? one : es(zero, one));
          oracle.emit(cell, E{sign == 1 ? half : mul(IOTA, half), 0, 0, 0});
        }
      }
      linear.constant =
          es(linear.constant, em(w, em(f.pub[s][0], f.pub[s][1 + j])));
      if (!fresh)
        linear.add_slot(f, s, 69 + j, es(zero, w));
    }
  for (int a = 5; a < 69; ++a) {
    E w = em(f.eta, power(f.rho, 9 + a - 5));
    linear.add_slot(f, f.k + 1, a, w);
    linear.add_slot(f, 0, a, es(zero, w));
    for (int s = 1; s <= f.k; ++s)
      linear.add_slot(f, s, a, es(zero, em(w, power(f.alpha, s))));
  }
  linear.constant = ea(linear.constant, one);
  oracle.emit(linear, E{half, 0, 0, 0});
  linear.constant = es(linear.constant, sc(one, 2));
  oracle.emit(linear, E{mul(IOTA, half), 0, 0, 0});
  check(oracle.values.size() == size_t(55 * V + 64 * (f.k + 2) + 2),
        "dense table size");
  const auto fast = fieldadj(f, eqweights(r));
  check(oracle.offset == fast.off, "dense public offset");
  check(oracle.row == fast.row, "dense affine row coefficients");
  E actual = fast.off;
  for (size_t j = 0; j < fast.row.size(); ++j)
    actual = ea(actual, sc(fast.row[j], cv(x, j)));
  check(actual == oracle.fold_table(), "non-Boolean table endpoint");
  return actual;
}

// Compare the original root residual to complete reconstructed physical
// openings. Public bits, inactive limbs and zero padding all contribute here.
inline void check_root(const PublicFront &f, const System &sys,
                       const Key &state, const Key &aux, const RV &x,
                       E endpoint) {
  SV residual;
  for (int s = 0; s < f.k + 2; ++s) {
    RV opening(219);
    for (int a = 0; a < 73; ++a)
      for (int t = 0; t < 4; ++t)
        for (int j = 0; j < 48; ++j) {
          int g = group(f.k, s, a, t);
          U digit = g >= 0  ? cv(x, 55 * g + j)
                    : a < 5 ? ((f.pub[s][a][t] >> j) & 1U)
                            : 0;
          put(opening, 48 * (4 * a + t) + j, digit);
        }
    RV image = state.apply(opening);
    for (size_t j = 0; j < image.size(); ++j)
      image[j] = rs(image[j], f.cs[s][j]);
    append(residual, image);
  }
  size_t groups = 64 * f.k + 544;
  RV digits(aux.cols);
  for (size_t i = 0; i < 7 * groups; ++i)
    put(digits, i, cv(x, 55 * (i / 7) + 48 + (i % 7)));
  RV image = aux.apply(digits);
  for (size_t j = 0; j < image.size(); ++j)
    image[j] = rs(image[j], f.caux[j]);
  append(residual, image);
  for (U v : es(endpoint, f.field.terminal))
    residual.push_back(v);
  residual.push_back(sub(cv(x, 55 * groups), 1));
  auto rows = sys.apply(x);
  check(rows.size() == residual.size(), "dense root row count");
  for (size_t j = 0; j < rows.size(); ++j)
    check(sub(rows[j], sys.y[j]) == residual[j], "dense full root residual");
}
} // namespace dense_handoff
