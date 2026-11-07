#include "sonata/epoch_session.hpp"
#include <iostream>
#include <stdexcept>

using namespace sonata;
namespace ep = sonata::epoch;

template<typename F> static void expect_rejection(F attempt) {
  try { attempt(); } catch (const std::runtime_error&) { return; }
  throw std::logic_error("expected rejection");
}

int main() {
  std::vector<E> point;
  for (const E v : {zero,one,E{9,8,0,1},one,zero,E{1,2,3,4}}) {
    auto expected=eqweights(point);
    ep::FactoredWeights weights(point);
    check(weights.size()==expected.size(),"split weight dimensions");
    for (size_t i=0;i<expected.size();++i)
      check(weights[i]==expected[i],"factor order or boundary challenge");
    bool threw=false;
    try { (void)weights[weights.size()]; } catch (const std::out_of_range&) { threw=true; }
    check(threw,"weight index outside support");
    point.push_back(v);
  }
  ep::Object state{}, src{};
  state.commitment.resize(4); src.commitment.resize(4); src.pub[0]=one;
  ep::Frozen frozen{{state,state},{zero,zero},{{{src},zero,one}},RV(4)};
  auto crs=std::make_shared<Master>(); crs->data.resize(44352);
  ep::AttemptMeter meter(2);
  size_t calls=0;
  // Reject one out-of-field limb and an all-zero candidate for rho,
  // then draw rho=eta=1 and zero Sumcheck points.
  auto coin=[&]() -> uint64_t {
    size_t i=calls++;
    if (i==0) return q;
    return i==5 || i==9 ? 1 : 0;
  };
  ep::FieldSession session(frozen,true,meter,coin,crs);
  ep::FieldSession aborted(frozen,true,meter,[]{return uint64_t(1);},crs);
  ep::FieldSession over_budget(frozen,true,meter,[]{return uint64_t(1);},crs);
  frozen.pending.back()=one; // Mutating the original cannot change the session.
  expect_rejection([&]{ ep::FieldSession bad(frozen,true,meter,coin,crs); });
  check(calls==0 && meter.remaining()==2,"invalid pre-batch prefix consumes no coin");
  expect_rejection([&]{ session.receive_polynomial(zero,zero); });
  auto batch=session.issue_batch();
  check(batch.first==one && batch.second==one && meter.remaining()==1,"post-freeze batch");
  expect_rejection([&]{ session.issue_batch(); });
  size_t rounds=0;
  for (size_t n=1;n<55*(64+2*272)+64*3+2;n*=2) ++rounds;
  for (size_t j=0;j<rounds;++j)
    check(session.receive_polynomial(zero,zero)==zero,"mock challenge schedule");
  // A consistent field recurrence produces a statement for the backend.
  // The all-zero opening below violates its fixed-one row.
  auto root=session.submit_terminal(zero);
  check(root.raw==3635520 && root.y.size()==264036,"physical handoff");
  check(root.apply(RV(root.raw/64))!=root.y,"no acceptance from mere field recurrence");
  expect_rejection([&]{ session.submit_terminal(zero); });
  aborted.issue_batch();
  aborted.abort();
  expect_rejection([&]{ aborted.receive_polynomial(zero,zero); });
  expect_rejection([&]{ over_budget.issue_batch(); });
  check(meter.remaining()==0,"abandoned attempt stays charged");
  std::cout<<"PASS: deep freeze, challenge order, handoff, abort and shared issue cap\n";
}
