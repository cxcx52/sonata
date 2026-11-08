#pragma once

#include "epoch_fixture.hpp"
#include "five_layer.hpp"

namespace sonata::example {

// The ordinary relation permits a nonzero fresh residual if the new public
// pending value accounts for it. The input fixture has exactly one such lane.
inline std::pair<System, RV> ordinary_transition(
    const std::shared_ptr<const Master>& crs, const WordSource& word) {
  constexpr int k = 1024;
  fixture::Fixture prover(k);
  PublicFront public_input{};
  public_input.k = k;
  Key state_key{crs, 4, 219};
  auto commit = [&](const fixture::State& state) {
    auto record = fixture::commit_state(state_key, state);
    public_input.pub.push_back(record.pub);
    public_input.cs.push_back(std::move(record.commitment));
  };
  std::mt19937_64 inputs(213);
  for (int j = 0; j <= k; ++j) {
    prover.states[j] = fixture::input_state(inputs, j != 0);
    if (j == 1)
      prover.states[j][1] = ea(prover.states[j][1], one);
    commit(prover.states[j]);
  }
  auto nonzero_coin = [&] {
    E value{};
    do { value = aggregation_coin(word); } while (value == zero);
    return value;
  };
  public_input.alpha = aggregation_coin(word);
  public_input.gamma = nonzero_coin();
  prover.gamma = public_input.gamma;
  prover.encode_fold(public_input.alpha);
  commit(prover.states.back());

  // The first source's x_1 was incremented by one while u=1. Its residual
  // is -1 and all remaining fresh residuals are zero.
  public_input.prev_pending = zero;
  public_input.pending = es(zero, public_input.gamma);
  const size_t groups = prover.scalarvalues.size();
  Key aux{crs, 4, (7 * groups + 63) / 64};
  RV opening(aux.cols);
  for (size_t g = 0; g < groups; ++g)
    for (size_t j = 0; j < 7; ++j)
      put(opening, 7*g+j, prover.bits[55*g+48+j]);
  public_input.caux = aux.apply(opening);

  public_input.rho = nonzero_coin();
  public_input.eta = nonzero_coin();
  prover.rho = public_input.rho;
  prover.eta = public_input.eta;
  prover.batch_equations();
  check(prover.invalid_sources == 1 && prover.pending == public_input.pending,
        "ordinary pending update");
  const E claim = ea(E{mul(groups, KAPPA), 0, 0, 0}, em(public_input.eta, prover.pending));
  auto received = fixture::prove_rounds(prover.table(), claim, true,
      [&](E g0, E g2) {
        public_input.field.g0.push_back(g0);
        public_input.field.g2.push_back(g2);
        E r = aggregation_coin(word);
        public_input.field.ch.push_back(r);
        return r;
      });
  public_input.field.terminal = received.terminal;
  check(public_input.field.ch.size() == 22 && groups == 66080, "ordinary profile");
  auto system = align(initialsystem(public_input, crs), layers.front());
  RV witness(layers.front().N / degree);
  for (size_t j = 0; j < prover.bits.size(); ++j)
    put(witness, j, prover.bits[j]);
  put(witness, prover.bits.size(), 1);
  std::cerr << "Ordinary root: k=1024 groups=" << groups
            << " sumcheck=22 nonzero_pending=1 rows=" << system.y.size() << '\n';
  return {std::move(system), std::move(witness)};
}

} // namespace sonata::example
