#pragma once
#include "sonata/epoch.hpp"
#include "sonata/epoch_envelope.hpp"

namespace sonata::epoch {

// A single sequential verifier owns this meter across retries and sessions.
class AttemptMeter {
  size_t left_;
public:
  explicit AttemptMeter(size_t limit) : left_(limit) {}
  AttemptMeter(const AttemptMeter&) = delete;
  AttemptMeter& operator=(const AttemptMeter&) = delete;
  void issue() {
    check(left_ != 0, "epoch attempt budget exhausted");
    --left_;
  }
  size_t remaining() const { return left_; }
};

// Starts with the auditor's retained log, after the auxiliary commitment.
// The terminal result is a backend statement; constructing it is not acceptance.
// A trusted outer verifier serializes the full attempt through the backend's
// verdict and holds the CRS immutable (including all aliases to the Master).
class FieldSession {
  enum class Phase { prepared, field, closed };
  Phase phase_ = Phase::prepared;
  const Frozen log_;
  const bool direct_;
  AttemptMeter& meter_;
  std::function<uint64_t()> word_;
  std::shared_ptr<const Master> crs_;
  size_t rounds_ = 0;
  E rho_{}, eta_{};
  Transcript transcript_;

  E draw(bool nonzero) {
    E v{};
    do {
      for (auto& a : v)
        do { a = word_() & mask; } while (a >= q);
    } while (nonzero && v == zero);
    return v;
  }
public:
  // word must supply independent uniform 64-bit words in the interactive model.
  FieldSession(const Frozen& frozen, bool direct, AttemptMeter& meter,
               std::function<uint64_t()> word, std::shared_ptr<const Master> crs)
    : log_(frozen), direct_(direct), meter_(meter), word_(std::move(word)),
      crs_(std::move(crs)) {
    auto l = layout(log_);
    check(l.groups <= 66080 && crs_ && crs_->data.size()==44352,
          "epoch outside fixed five-layer profile");
    if (direct_)
      for (size_t t=1; t<log_.pending.size(); ++t)
        check(log_.pending[t]==log_.pending[t-1],"direct pending continuity");
    for (size_t t=1; t<l.table_size(); t*=2) ++rounds_;
  }
  FieldSession(const FieldSession&) = delete;
  FieldSession& operator=(const FieldSession&) = delete;

  std::pair<E,E> issue_batch() {
    check(phase_==Phase::prepared,"batch coins already issued or session closed");
    phase_=Phase::closed;
    meter_.issue();
    rho_=draw(true);
    eta_=draw(true);
    phase_=Phase::field;
    return {rho_,eta_};
  }

  E receive_polynomial(E g0,E g2) {
    check(phase_==Phase::field,"no active field round");
    phase_=Phase::closed;
    check(canonical(g0) && canonical(g2) && transcript_.ch.size()<rounds_,
          "malformed or extra sumcheck round");
    transcript_.g0.push_back(g0);
    transcript_.g2.push_back(g2);
    E coin=draw(false);
    transcript_.ch.push_back(coin);
    phase_=Phase::field;
    return coin;
  }

  System submit_terminal(E endpoint) {
    check(phase_==Phase::field,"terminal outside live field phase");
    phase_=Phase::closed;
    check(transcript_.ch.size()==rounds_,"incomplete sumcheck");
    transcript_.terminal=endpoint;
    const Config cfg{21,2705,64,3634401,3635520,3105761,795074815,1174,22};
    auto active=initialsystem(log_,rho_,eta_,transcript_,crs_,7228,direct_);
    return align(embed_prefix(std::move(active),cfg,264036),cfg);
  }
  void abort() { phase_=Phase::closed; }
};
} // namespace sonata::epoch
