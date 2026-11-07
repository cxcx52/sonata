#pragma once
#include "epoch_fixture.hpp"
#include "sonata/epoch_session.hpp"

namespace sonata::example {
struct WindowInput {
  System system;
  RV witness;
};

// Starts with an honest history rooted at a trusted synthetic anchor. Only
// the honest prover uses the opening. FieldSession receives the public log.
inline WindowInput make_window_input(std::shared_ptr<const Master> master,
                                     std::mt19937_64& field_rng) {
  check(master && master->data.size() == 44352, "expected the full-profile CRS");
  auto [log, witness] = fixture::make_window_opening(
      master, fixture::target_window_arities(), zero);
  const auto layout = epoch::layout(log);
  check(log.steps.size() == 8 && layout.blocks.size() - log.states.size() == 994 &&
            layout.groups == 66064 && log.pending.back() == zero,
        "expected eight-step zero-pending window");
  fixture::commit_window_aux(log, witness, Key{master, 4, 7228}, layout.groups);

  epoch::AttemptMeter attempts(1);
  constexpr bool direct = false; // Retain the pending equation at every step.
  epoch::FieldSession verifier(log, direct, attempts,
                               [&field_rng] { return field_rng(); }, master);
  const auto [rho, eta] = verifier.issue_batch();
  auto table = fixture::forward(log, rho, eta, witness, direct);
  const E claim{mul(layout.groups, KAPPA), 0, 0, 0};
  auto transcript = fixture::prove_rounds(std::move(table), claim, true,
      [&verifier](E g0, E g2) { return verifier.receive_polynomial(g0, g2); });
  System system = verifier.submit_terminal(transcript.terminal);
  check(transcript.ch.size() == 22 && system.raw == 3635520 &&
            system.y.size() == 264036 && attempts.remaining() == 0,
        "unexpected frozen-window handoff dimensions");
  witness.resize(system.raw / 64);
  return {std::move(system), std::move(witness)};
}
} // namespace sonata::example
