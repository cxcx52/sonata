#pragma once

#include "epoch_fixture.hpp"
#include "five_layer.hpp"

namespace sonata::example {

// The ordinary relation permits a nonzero fresh residual if the new public
// pending value accounts for it. The input fixture has exactly one such lane.
struct OrdinaryInputs {
  fixture::Fixture prover;
  PublicFront pub;
};

inline OrdinaryInputs prepare_ordinary_inputs(int k, const std::shared_ptr<const Master>& crs) {
  OrdinaryInputs inputs{fixture::Fixture(k), {}};
  inputs.pub.k = k;
  inputs.pub.prev_pending = zero;
  std::mt19937_64 rng(213);
  Key state_key{crs, 4, 219};
  for (int j = 0; j <= k; ++j) {
    auto &state = inputs.prover.states[j];
    state = fixture::input_state(rng, j != 0);
    if (j == 1) state[1] = ea(state[1], one);
    auto record = fixture::commit_state(state_key, state);
    inputs.pub.pub.push_back(record.pub);
    inputs.pub.cs.push_back(std::move(record.commitment));
  }
  return inputs;
}

inline std::pair<System, RV> finish_ordinary_transition(
    OrdinaryInputs &input, const std::shared_ptr<const Master>& crs, const WordSource& word,
    const Config& root, double *coin_seconds = nullptr) {
  auto &public_input = input.pub;
  auto &prover = input.prover;
  const int k = public_input.k;
  auto nonzero_coin = [&] {
    E value{};
    do { value = aggregation_coin(word); } while (value == zero);
    return value;
  };
  public_input.alpha = account(coin_seconds, [&]{return aggregation_coin(word);});
  public_input.gamma = account(coin_seconds, nonzero_coin);
  prover.gamma = public_input.gamma;
  prover.encode_fold(public_input.alpha);
  auto record = fixture::commit_state(Key{crs, 4, 219}, prover.states.back());
  public_input.pub.push_back(record.pub);
  public_input.cs.push_back(std::move(record.commitment));

  // The first source's x_1 was incremented by one while u=1. Its residual
  // is -1 and all remaining fresh residuals are zero.
  public_input.pending = es(zero, public_input.gamma);
  const size_t groups = prover.scalarvalues.size();
  Key aux{crs, 4, (7 * groups + 63) / 64};
  RV opening(aux.cols);
  for (size_t g = 0; g < groups; ++g)
    for (size_t j = 0; j < 7; ++j)
      put(opening, 7*g+j, prover.bits[55*g+48+j]);
  public_input.caux = aux.apply(opening);

  public_input.rho = account(coin_seconds, nonzero_coin);
  public_input.eta = account(coin_seconds, nonzero_coin);
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
        E r = account(coin_seconds, [&]{return aggregation_coin(word);});
        public_input.field.ch.push_back(r);
        return r;
      });
  public_input.field.terminal = received.terminal;
  check(groups == size_t(64*k + 544) && root.raw == 55*groups+1, "ordinary profile");
  auto system = align(initialsystem(public_input, crs), root);
  RV witness(root.N / degree);
  for (size_t j = 0; j < prover.bits.size(); ++j)
    put(witness, j, prover.bits[j]);
  put(witness, prover.bits.size(), 1);
  std::cerr << "Ordinary root: k=" << k << " groups=" << groups
            << " sumcheck=" << public_input.field.ch.size() << " nonzero_pending=1 rows=" << system.y.size() << '\n';
  return {std::move(system), std::move(witness)};
}

inline std::pair<System, RV> ordinary_transition(
    const std::shared_ptr<const Master>& crs, const WordSource& word) {
  auto input = prepare_ordinary_inputs(1024, crs);
  return finish_ordinary_transition(input, crs, word, layers.front());
}

} // namespace sonata::example
