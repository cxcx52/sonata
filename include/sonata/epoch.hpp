#pragma once
#include "sonata/ring.hpp"
#include "sonata/epoch_weights.hpp"
#include <limits>

namespace sonata::epoch {
struct Object {
  std::array<E, 5> pub;
  RV commitment;
};
struct Step {
  std::vector<Object> sources;
  E alpha, gamma;
};
struct Frozen {
  std::vector<Object> states;
  std::vector<E> pending;
  std::vector<Step> steps;
  RV aux_commitment;
};
struct Block {
  size_t first;
  bool fresh;
};
struct Layout {
  std::vector<Block> blocks;
  size_t groups = 0;
  size_t table_size() const { return 55 * groups + 64 * blocks.size() + 2; }
  size_t raw() const { return 55 * groups + 1; }
};
inline bool canonical(const E &x) {
  return std::all_of(x.begin(), x.end(), [](U v) { return v < q; });
}
inline void check_commitment(const RV &c) {
  check(c.size() == 4, "epoch commitment height");
  for (const auto &r : c)
    check(std::all_of(r.begin(), r.end(), [](U a) { return a < q; }),
          "noncanonical commitment");
}
inline Layout layout(const Frozen &f) {
  check(!f.steps.empty() && f.states.size() == f.steps.size() + 1 &&
            f.pending.size() == f.states.size(),
        "epoch state count");
  Layout l;
  auto add_object = [&](const Object &o, bool fresh) {
    check_commitment(o.commitment);
    for (const auto &e : o.pub)
      check(canonical(e), "noncanonical public coordinate");
    if (fresh) {
      check(o.pub[0] == one, "fresh public u");
      for (size_t i = 1; i < 5; i++)
        check(o.pub[i][1] == 0 && o.pub[i][2] == 0 && o.pub[i][3] == 0,
              "fresh public field");
    }
    check(l.groups < std::numeric_limits<size_t>::max() / 256 - 272,
          "epoch size overflow");
    l.blocks.push_back({l.groups, fresh});
    l.groups += fresh ? 64 : 272;
  };
  for (const auto &s : f.states)
    add_object(s, false);
  for (const auto &p : f.pending)
    check(canonical(p), "noncanonical pending");
  for (size_t i = 0; i < f.steps.size(); ++i) {
    const auto &t = f.steps[i];
    check(!t.sources.empty() && canonical(t.alpha) && canonical(t.gamma) &&
              t.gamma != zero,
          "epoch step fields");
    auto folded = f.states[i].pub;
    E power = one;
    for (const auto &src : t.sources) {
      add_object(src, true);
      power = em(power, t.alpha);
      for (size_t j = 0; j < 5; ++j)
        folded[j] = ea(folded[j], em(power, src.pub[j]));
    }
    check(folded == f.states[i + 1].pub, "epoch public state link");
  }
  check_commitment(f.aux_commitment);
  return l;
}

// Logical variables are the decoded private slots 5..72. Product rows use
// consecutive pairs among 5..68. The final affine form also refers to errors.
struct Term {
  size_t object, slot;
  E weight;
};
struct Batch {
  std::vector<std::array<E, 4>> products;
  std::vector<Term> linear;
  E constant{};
};
inline Batch batch(const Frozen &f, const Layout &l, E rho, E eta,
                   bool drop_pending_equations = false) {
  check(canonical(rho) && canonical(eta) && rho != zero && eta != zero,
        "epoch batch coins");
  Batch b;
  b.products.resize(l.blocks.size());
  E power = eta;
  auto next = [&]() {
    power = em(power, rho);
    return power;
  };
  auto quad = [&](size_t i, size_t j, const Object &o, E weight) {
    b.products[i][j] = weight;
    b.constant = es(b.constant, em(weight, em(o.pub[0], o.pub[j + 1])));
  };
  // Equation order: all state relations, then for each step its fresh
  // relations, its private folds and its pending update (new minus old).
  for (size_t i = 0; i < f.states.size(); ++i)
    for (size_t j = 0; j < 4; ++j) {
      E w = next();
      quad(i, j, f.states[i], w);
      b.linear.push_back({i, 69 + j, es(zero, w)});
    }
  size_t first = f.states.size();
  for (size_t i = 0; i < f.steps.size(); ++i) {
    const auto &t = f.steps[i];
    for (size_t k = 0; k < t.sources.size(); ++k)
      for (auto &w : b.products[first + k])
        w = next();
    for (size_t slot = 5; slot < 69; ++slot) {
      E w = next(), a = w;
      b.linear.push_back({i + 1, slot, w});
      b.linear.push_back({i, slot, es(zero, w)});
      for (size_t k = 0; k < t.sources.size(); ++k) {
        a = em(a, t.alpha);
        b.linear.push_back({first + k, slot, es(zero, a)});
      }
    }
    if (drop_pending_equations)
      check(f.pending[i] == f.pending[i + 1],
            "direct settlement pending continuity");
    E update = drop_pending_equations ? zero : next(), g = update;
    b.constant = ea(b.constant, em(update, es(f.pending[i + 1], f.pending[i])));
    for (size_t k = 0; k < t.sources.size(); ++k)
      for (size_t j = 0; j < 4; ++j) {
        g = em(g, t.gamma);
        quad(first + k, j, t.sources[k], es(b.products[first + k][j], g));
      }
    first += t.sources.size();
  }
  return b;
}

struct Endpoint {
  std::vector<E> row;
  E offset{};
};
inline Endpoint endpoint(const Layout &l, const Batch &b,
                         const std::vector<E> &r) {
#ifdef SONATA_EPOCH_DENSE_WEIGHTS
  const auto ew = eqweights(r);
#else
  const FactoredWeights ew(r);
#endif
  check(ew.size() >= l.table_size(), "epoch endpoint weights");
  Endpoint out;
  out.row.resize(55 * l.groups);
  size_t pos = 0;
  for (size_t g = 0; g < l.groups; ++g) {
    E prefix{};
    for (size_t i = 0; i < 55; ++i) {
      E rt{SQRT_D[i][0], SQRT_D[i][1], SQRT_D[i][2], SQRT_D[i][3]};
      E a = em(ew[pos++], rt);
      out.offset = ea(out.offset, sc(a, BETA[i]));
      if (i >= 6 && i < 48) {
        out.row[55 * g + i] = ea(out.row[55 * g + i], ea(a, prefix));
        prefix = ea(prefix, sc(a, high_inv[i - 6]));
      } else
        for (size_t j = i; j < 55; ++j)
          if (LDL_L[j][i])
            out.row[55 * g + j] = ea(out.row[55 * g + j], sc(a, LDL_L[j][i]));
    }
    for (size_t j = 0; j < 5; ++j)
      out.row[55 * g + 48 + j] =
          ea(out.row[55 * g + 48 + j], sc(prefix, 2 * HWEIGHTS[j]));
    out.row[55 * g + 53] = es(out.row[55 * g + 53], prefix);
  }
  // Accumulate coefficients in decoded K slots first, then expand each
  // coordinate once into its 48 value limbs.
  std::vector<std::array<E, 68>> slot_weights(l.blocks.size());
  auto pair_weights = [&]() {
    E u = ew[pos++];
    E v = sc(ew[pos++], IOTA);
    return std::make_pair(sc(ea(u, v), iv2), sc(es(u, v), iv2));
  };
  for (size_t obj = 0; obj < l.blocks.size(); ++obj)
    for (size_t j = 0; j < 32; ++j) {
      auto [a, c] = pair_weights();
      slot_weights[obj][2 * j] = em(a, b.products[obj][j / 8]);
      slot_weights[obj][2 * j + 1] = c;
    }
  auto [a, c] = pair_weights();
  out.offset = ea(out.offset, ea(c, em(a, b.constant)));
  for (const auto &term : b.linear) {
    auto &w = slot_weights[term.object][term.slot - 5];
    w = ea(w, em(a, term.weight));
  }
  for (size_t obj = 0; obj < l.blocks.size(); ++obj) {
    auto block = l.blocks[obj];
    size_t dims = block.fresh ? 1 : 4;
    for (size_t j = 0; j < (block.fresh ? 64u : 68u); ++j)
      for (size_t t = 0; t < dims; ++t) {
        E base{};
        base[t] = 1;
        E w = em(slot_weights[obj][j], base);
        size_t start = 55 * (block.first + dims * j + t);
        for (size_t bit = 0; bit < 48; ++bit) {
          out.row[start + bit] = ea(out.row[start + bit], w);
          w = ea(w, w);
        }
      }
  }
  return out;
}

// Map each active scalar to its original compressed state's physical bit
// position. Public/reserved slots do not occur in the free vector.
template <class F> inline void each_value(const Block &b, F visit) {
  size_t dims = b.fresh ? 1 : 4;
  for (size_t j = 0; j < (b.fresh ? 64u : 68u); ++j)
    for (size_t t = 0; t < dims; ++t)
      for (size_t bit = 0; bit < 48; ++bit)
        visit(55 * (b.first + dims * j + t) + bit,
              48 * (4 * (j + 5) + t) + bit);
}
inline System initialsystem(const Frozen &f, E rho, E eta, const Transcript &tr,
                            std::shared_ptr<const Master> master,
                            size_t aux_columns = 0,
                            bool drop_pending_equations = false) {
  auto l = std::make_shared<const Layout>(layout(f));
  size_t ell = 0;
  for (size_t n = l->table_size() - 1; n; n >>= 1)
    ++ell;
  check(tr.ch.size() == ell && tr.g0.size() == ell && tr.g2.size() == ell,
        "epoch sumcheck shape");
  for (const auto *v : {&tr.ch, &tr.g0, &tr.g2})
    for (const auto &e : *v)
      check(canonical(e), "noncanonical sumcheck field");
  check(canonical(tr.terminal) &&
            verify(tr, E{mul(l->groups % q, KAPPA), 0, 0, 0}, tr.terminal),
        "epoch sumcheck recurrence");
  auto ep = std::make_shared<const Endpoint>(
      endpoint(*l, batch(f, *l, rho, eta, drop_pending_equations), tr.ch));
  size_t min_aux = (7 * l->groups + 63) / 64;
  if (!aux_columns)
    aux_columns = min_aux;
  check(master && aux_columns >= min_aux &&
            std::max(aux_columns, size_t(219)) <= master->data.size() / 4,
        "epoch CRS prefix width");
  Key st{master, 4, 219}, aux{master, 4, aux_columns};
  System sys;
  sys.raw = l->raw();
  auto bind = [&](const Object &obj) {
    RV fixed(219);
    for (size_t j = 0; j < 5; ++j)
      for (size_t t = 0; t < 4; ++t)
        for (size_t bit = 0; bit < 48; ++bit)
          put(fixed, 48 * (4 * j + t) + bit, (obj.pub[j][t] >> bit) & 1U);
    auto offset = st.apply(fixed);
    for (size_t j = 0; j < 4; ++j)
      offset[j] = rs(obj.commitment[j], offset[j]);
    append(sys.y, offset);
  };
  for (const auto &obj : f.states)
    bind(obj);
  for (const auto &step : f.steps)
    for (const auto &obj : step.sources)
      bind(obj);
  append(sys.y, f.aux_commitment);
  for (U z : es(tr.terminal, ep->offset))
    sys.y.push_back(z);
  sys.y.push_back(1);
  sys.apply = [l, ep, st, aux](const RV &x) {
    SV y;
    for (const auto &b : l->blocks) {
      RV packed(st.cols);
      each_value(b, [&](size_t j, size_t p) { put(packed, p, cv(x, j)); });
      append(y, st.apply(packed));
    }
    RV packed(aux.cols);
    for (size_t j = 0; j < 7 * l->groups; ++j)
      put(packed, j, cv(x, 55 * (j / 7) + 48 + j % 7));
    append(y, aux.apply(packed));
    E v{};
    for (size_t j = 0; j < ep->row.size(); ++j)
      v = ea(v, sc(ep->row[j], cv(x, j)));
    for (U z : v)
      y.push_back(z);
    y.push_back(cv(x, l->raw() - 1));
    return y;
  };
  const size_t rows = sys.y.size();
  sys.adj = [l, ep, st, aux, rows](const SV &y) {
    check(y.size() == rows, "epoch adjoint dimensions");
    RV out((l->raw() + 63) / 64);
    size_t p = 0;
    auto transpose = [&](const Key &key) {
      auto v = key.transpose(chunks(y, p, 4, true));
      p += 256;
      for (auto &cell : v)
        cell = dual(cell);
      return v;
    };
    for (const auto &b : l->blocks) {
      auto v = transpose(st);
      each_value(b, [&](size_t j, size_t phys) { put(out, j, cv(v, phys)); });
    }
    auto v = transpose(aux);
    for (size_t j = 0; j < 7 * l->groups; ++j)
      put(out, 55 * (j / 7) + 48 + j % 7, cv(v, j));
    for (size_t j = 0; j < ep->row.size(); ++j) {
      U a = cv(out, j);
      for (size_t t = 0; t < 4; ++t)
        a = add(a, mul(y[p + t], ep->row[j][t]));
      put(out, j, a);
    }
    put(out, l->raw() - 1, y[p + 4]);
    for (auto &cell : out)
      cell = dual(cell);
    return out;
  };
  return sys;
}
} // namespace sonata::epoch
