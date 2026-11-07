#pragma once
#include "sonata/epoch_session.hpp"

inline void check_frozen_honest_handoff(const sonata::epoch::Frozen& f,
    sonata::E rho,sonata::E eta,const sonata::Transcript& tr,sonata::RV w,
    const std::shared_ptr<sonata::Master>& test_master) {
  using namespace sonata;
  auto crs=std::make_shared<Master>(*test_master);
  crs->data.resize(44352); // Prefix and every used opening are unchanged.
  SV tape;
  for(const auto& v:{rho,eta}) tape.insert(tape.end(),v.begin(),v.end());
  for(const auto& v:tr.ch) tape.insert(tape.end(),v.begin(),v.end());
  size_t at=0;
  epoch::AttemptMeter attempts(1);
  auto original=f;
  epoch::FieldSession verifier(original,true,attempts,[&]{return tape.at(at++);},crs);
  auto coins=verifier.issue_batch();
  check(coins==std::make_pair(rho,eta),"fixture batch challenges");
  original.pending.back()=ea(original.pending.back(),one);
  for(size_t j=0;j<tr.ch.size();++j)
    check(verifier.receive_polynomial(tr.g0[j],tr.g2[j])==tr.ch[j],"field round order");
  auto physical=verifier.submit_terminal(tr.terminal);
  w.resize(physical.raw/64);
  check(at==tape.size() && attempts.remaining()==0 && physical.apply(w)==physical.y,
        "original frozen root after external mutation");
  std::cout<<"PASS: honest frozen log, every field message and original physical root\n";
}
