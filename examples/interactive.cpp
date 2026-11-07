// Live, small-profile strict checkpoint. The prover owns fixture data and
// openings. Only the public receive methods below issue verifier challenges.
#include "challenge_fixture.hpp"
#include "five_layer.hpp"
#include "epoch_fixture.hpp"
#include "sonata/epoch_session.hpp"
#include <fstream>

namespace {
using namespace sonata;
using namespace sonata::example;
namespace ep = sonata::epoch;
namespace fx = sonata::fixture;

// Two steps with arities 1 and 2. The second S follows child_layout(289,47377,3,64).
constexpr Config root_cfg{3,289,64,55441,55488,47377,256ULL*47377-1,1120,16};
constexpr Config last_cfg{3,385,0,73856,73920,56413322,256ULL*56413322,1152,16};
enum class Fault { none, commitment, residual, pending, continuity, response };

// Only a failure thrown from an explicitly marked public receiver can
// satisfy a negative test. Every exception ends the current interaction.
class VerifierRejection : public std::runtime_error {
public:
  const char* stage;
  VerifierRejection(const char* where, const char* why)
    : std::runtime_error(why), stage(where) {}
};
template <class F> decltype(auto) receive(const char* stage, F&& fn) {
  try { return fn(); }
  catch (const std::runtime_error& e) { throw VerifierRejection(stage,e.what()); }
}
struct Issued {
  size_t alpha=0, gamma=0, batches=0, rounds=0;
};

class VerifierCoins {
  std::ifstream entropy_{"/dev/urandom", std::ios::binary};
public:
  U word() {
    U result;
    if (!entropy_.read(reinterpret_cast<char*>(&result), sizeof(result)))
      throw std::runtime_error("verifier entropy unavailable");
    return result;
  }
  E field(bool nonzero) {
    E value;
    do { value = aggregation_coin([this] { return word(); }); }
    while (nonzero && value == zero);
    return value;
  }
};

// Collector starts at a trusted anchor and retains the public log by value.
class Collector {
  VerifierCoins& coins_;
  Issued& issued_;
  ep::Frozen log_;
public:
  Collector(VerifierCoins& coins, Issued& issued, ep::Object anchor) : coins_(coins), issued_(issued) {
    log_.states.push_back(std::move(anchor));
    log_.pending.push_back(zero);
  }
  std::pair<E,E> begin(const ep::Object& old, E pending, std::vector<ep::Object> sources) {
    check(log_.steps.size()+1 == log_.states.size(), "unfinished transition");
    const auto& expected = log_.states.back();
    check(old.pub == expected.pub && old.commitment == expected.commitment &&
              pending == log_.pending.back(), "state continuity");
    check(!sources.empty(), "empty step");
    for (const auto& source : sources) ep::check_commitment(source.commitment);
    log_.steps.push_back({std::move(sources), {}, {}});
    auto& step = log_.steps.back();
    step.alpha = coins_.field(false);
    ++issued_.alpha;
    step.gamma = coins_.field(true);
    ++issued_.gamma;
    return {step.alpha, step.gamma};
  }
  void end(ep::Object output, E pending) {
    check(log_.steps.size() == log_.states.size(), "no open transition");
    ep::check_commitment(output.commitment);
    log_.states.push_back(std::move(output));
    log_.pending.push_back(pending);
  }
  ep::Frozen freeze(RV auxiliary) && {
    check(log_.states.size() == log_.steps.size()+1, "unfinished window");
    ep::check_commitment(auxiliary);
    log_.aux_commitment = std::move(auxiliary);
    ep::layout(log_);
    return std::move(log_);
  }
};

// The public receiver binds each message before revealing the next coin.
// A spent retry cannot be revisited after the receiver advances.
class BackendVerifier {
  enum Phase { seeds, image, auxiliary, short_coins, answered, closed };
  Phase phase_ = seeds;
  Config cfg_;
  System system_;
  LayerKeys keys_;
  Key pivot_key_;
  VerifierCoins& coins_;
  Bytes seed_;
  RV u1_, u2_, c_;
  TerminalResponse terminal_;
  Aggregated ag_;
  unsigned seed_count_ = 0, zeta_count_ = 0, challenge_count_ = 0;
public:
  BackendVerifier(Config cfg, System sys, const std::shared_ptr<const Master>& crs,
                  VerifierCoins& coins, RV tag, RV other_blocks)
      : cfg_(cfg), system_(std::move(sys)), keys_(cfg,crs),
        pivot_key_{crs,2,256}, coins_(coins) {
    check_canonical_rings(tag);
    check_canonical_rings(other_blocks);
    if (cfg.rho) {
      check(tag.size()==6 && other_blocks.empty(), "pre-projection commitment");
      u1_=std::move(tag);
    } else {
      check(tag.size()==2 && other_blocks.size()==16*(cfg.s-1), "pivot precommitment");
      terminal_.pivot=std::move(tag);
      terminal_.tothers=std::move(other_blocks);
    }
  }
  Bytes issue_seed() {
    check((phase_==seeds || phase_==image) && seed_count_ < attempt_limit, "seed retry");
    ++seed_count_;
    seed_=projection_seed(cfg_, [this]{return coins_.word();});
    phase_=image;
    return seed_;
  }
  E receive_image(const SV& p) {
    check(phase_==image && std::all_of(p.begin(),p.end(),[](U v){return v<q;}) &&
          proj_norm(p,cfg_.S), "projection image");
    phase_=auxiliary;
    E zeta=coins_.field(false);
    ++zeta_count_;
    ag_=aggregate(system_,expand_projection(cfg_,seed_),p,zeta,one);
    return zeta;
  }
  void receive_aux(RV msg) {
    check(phase_==auxiliary,"auxiliary order");
    check_canonical_rings(msg);
    if (cfg_.rho) {
      check(msg.size()==6,"child auxiliary commitment");
      u2_=std::move(msg);
    } else {
      const size_t triangle=cfg_.s*(cfg_.s+1)/2;
      check(msg.size()==4*triangle,"pivot h dimension");
      for (size_t e=0; e<4; ++e) check(msg[e*triangle][0]==0,"pivot h omitted coefficient");
      terminal_.h=std::move(msg);
    }
    phase_=short_coins;
  }
  RV issue_short() {
    check((phase_==short_coins || phase_==answered) &&
          challenge_count_ < attempt_limit,"short retry");
    ++challenge_count_;
    c_=draw_challenges(cfg_.s,!cfg_.rho);
    phase_=answered;
    return c_;
  }
  System commit_child() {
    check(cfg_.rho && phase_==answered,"child submission order");
    phase_=closed;
    return align(childsystem({cfg_,keys_.A,keys_.B,keys_.D,u1_,u2_,c_,std::move(ag_),
                              keys_.B.cols,keys_.D.cols,keys_.ell}),last_cfg);
  }
  void finish(RV z) {
    check(!cfg_.rho && phase_==answered,"terminal submission order");
    phase_=closed;
    check_canonical_rings(z);
    terminal_.z=std::move(z);
    verifyterminal(cfg_,keys_.A,pivot_key_,ag_,terminal_,c_);
  }
  void trace() const {
    std::cout << " seed=" << seed_count_ << " zeta=" << zeta_count_
              << " short=" << challenge_count_ << '\n';
  }
};

void backend(System prover, System verifier, RV w, const std::shared_ptr<const Master>& crs,
             VerifierCoins& coins, Fault fault) {
  for (const Config cfg : {root_cfg,last_cfg}) {
    check(w.size()*64==cfg.N && energy(w)<=cfg.S && prover.apply(w)==prover.y,
          "prover opening construction");
    LayerKeys keys(cfg,crs);
    RV tags;
    for (size_t j=0; j<cfg.s; ++j) {
      auto block=keys.A.apply(slice(w,j*cfg.n,cfg.n));
      tags.insert(tags.end(),block.begin(),block.end());
    }
    RV tag,bt,others;
    if (cfg.rho) {
      bt=digits(tags,cfg.rho,keys.ell);
      tag=keys.B.apply(bt);
    } else {
      tag=Key{crs,2,256}.apply(digits(slice(tags,0,16),8,16));
      others.assign(tags.begin()+16,tags.end());
    }
    BackendVerifier live(cfg,std::move(verifier),crs,coins,tag,std::move(others));
    Aggregated ag;
    for (;;) {
      auto P=expand_projection(cfg,live.issue_seed());
      auto p=project(P,w);
      if (proj_norm(p,cfg.S)) {
        E zeta=live.receive_image(p);
        ag=aggregate(prover,P,p,zeta,one);
        break;
      }
    }
    auto h=symmetric(ag,w,cfg.s,cfg.n);
    RV bh,u2;
    if (cfg.rho) {
      bh=digits(h,cfg.rho,keys.ell);
      u2=keys.D.apply(bh);
      live.receive_aux(u2);
    } else {
      const size_t triangle=cfg.s*(cfg.s+1)/2;
      for (size_t e=0; e<4; ++e) h[e*triangle][0]=0;
      live.receive_aux(std::move(h));
    }
    RV c,z;
    do {
      c=live.issue_short();
      z=response(w,c,cfg.n);
    } while (energy(z)>cfg.G);
    std::cout << " layer=" << (cfg.rho ? 1 : 2);
    live.trace();
    if (cfg.rho) {
      verifier=live.commit_child();
      prover=align(childsystem({cfg,keys.A,keys.B,keys.D,tag,u2,c,std::move(ag),
                                keys.B.cols,keys.D.cols,keys.ell}),last_cfg);
      w=child_opening(z,cfg.rho,bt,bh);
      w.resize(last_cfg.N/64);
    } else {
      if (fault==Fault::response) {
        // Preserve the response norm so rejection reaches the pivot equation.
        bool changed=false;
        for (auto& ring : z) {
          auto it=std::find_if(ring.begin(),ring.end(),[](U v){return v!=0;});
          if (it!=ring.end()) { *it=q-*it; changed=true; break; }
        }
        if (!changed) z[0][0]=1; // A unit vector also passes this profile's cap.
      }
      receive("pivot",[&]{live.finish(std::move(z));});
    }
  }
}

fx::State folded(const fx::State& old, const std::vector<fx::State>& sources, E alpha) {
  fx::State out=old;
  E power=one;
  for (const auto& source : sources) {
    power=em(power,alpha);
    for (size_t i=0; i<69; ++i) out[i]=ea(out[i],em(power,source[i]));
  }
  for (size_t lane=0; lane<4; ++lane) {
    E residual=es(zero,em(out[0],out[1+lane]));
    for (size_t j=5+16*lane; j<5+16*(lane+1); j+=2)
      residual=ea(residual,em(out[j],out[j+1]));
    out[69+lane]=residual;
  }
  return out;
}

void run_case(const std::shared_ptr<const Master>& crs, Fault fault) {
  VerifierCoins coins;
  std::mt19937_64 inputs(213); // Public test data only.
  Key st{crs,4,219};
  std::vector<fx::State> states{fx::input_state(inputs,false)}, sources;
  auto current=fx::commit_state(st,states.back());
  Issued issued;
  Collector collector(coins,issued,current);
  ep::Frozen honest;
  E pending=zero;
  honest.states.push_back(current);
  honest.pending.push_back(pending);
  for (size_t step=0; step<2; ++step) {
    std::vector<fx::State> batch;
    std::vector<ep::Object> public_sources;
    for (size_t i=0; i<=step; ++i) {
      auto source=fx::input_state(inputs,true);
      // Change x before commitment, keeping its private factors. Since u=1,
      // the first fresh residual becomes exactly -1.
      if (fault==Fault::residual && step==0) source[1]=ea(source[1],one);
      public_sources.push_back(fx::commit_state(st,source));
      batch.push_back(source);
      sources.push_back(source);
    }
    if (fault==Fault::continuity && step==1)
      current.commitment[0][0]=add(current.commitment[0][0],1);
    auto [alpha,gamma]=receive("collection",[&]{return collector.begin(current,pending,public_sources);});
    states.push_back(folded(states.back(),batch,alpha));
    current=fx::commit_state(st,states.back());
    honest.steps.push_back({public_sources,alpha,gamma});
    honest.states.push_back(current);
    honest.pending.push_back(zero);
    // A fabricated intermediate pending increment is used as the matching
    // old pending on the next step; both ends of the window are still zero.
    pending=(fault==Fault::pending && step==0) ? one : zero;
    collector.end(current,pending);
  }
  states.insert(states.end(),sources.begin(),sources.end());
  honest.aux_commitment.resize(4);
  auto w=fx::window_witness(honest,states);
  auto layout=ep::layout(honest);
  const size_t auxcols=(7*layout.groups+63)/64;
  fx::commit_window_aux(honest,w,Key{crs,4,auxcols},layout.groups);
  auto tag=honest.aux_commitment;
  if (fault==Fault::commitment) tag[0][0]=add(tag[0][0],1);
  auto public_log=std::move(collector).freeze(std::move(tag));
  ep::AttemptMeter meter(1);
  ep::FieldSession live(public_log,false,meter,[&coins]{return coins.word();},crs,auxcols);
  auto [rho,eta]=live.issue_batch();
  ++issued.batches;
  auto table=fx::forward(honest,rho,eta,w,false);
  // A prover with an invalid source can still send a complete transcript
  // for its actual table sum. The receiver enforces its own initial claim.
  E table_sum{};
  for (const auto& entry : table) table_sum=ea(table_sum,sq(entry));
  auto tr=fx::prove_rounds(std::move(table),table_sum,true,
      [&](E g0,E g2){
        auto r=live.receive_polynomial(g0,g2);
        ++issued.rounds;
        return r;
      });
  check(meter.remaining()==0 && issued.alpha==2 && issued.gamma==2 && issued.batches==1 &&
        issued.rounds==16 && issued.rounds==tr.ch.size(),"field challenge count");
  std::cout << " online_alpha=" << issued.alpha << " online_gamma=" << issued.gamma
            << " batch_pairs=" << issued.batches << " sumcheck=" << issued.rounds << '\n';
  auto verifier_root=align(receive("field",[&]{return live.submit_active_terminal(tr.terminal);}),root_cfg);
  auto prover_root=align(ep::initialsystem(honest,rho,eta,tr,crs,auxcols),root_cfg);
  w.resize(root_cfg.N/64);
  backend(std::move(prover_root),std::move(verifier_root),std::move(w),crs,coins,fault);
  check(public_log.pending.front()==zero && public_log.pending.back()==zero,
        "expected zero-pending decision");
}

} // namespace

int main(int argc,char** argv) {
  if (argc!=2 || std::string_view(argv[1])!="--small-live" ||
      !std::getenv("SONATA_E2E_CHALLENGE_PIPE")) {
    std::cerr << "Usage: python3 examples/run_interactive.py ./build/sonata_interactive\n";
    return 2;
  }
  try {
    auto master=std::make_shared<sonata::Master>();
    std::mt19937_64 setup(818); // Common setup, independent of the verifier coins.
    master->data.resize(44352);
    for (auto& r : master->data) r=sonata::fixture::randR(setup);
    struct Case {Fault fault; const char* name; const char* stage; const char* rejection;};
    for (const Case test : {Case{Fault::none,"honest","",""},
                             {Fault::commitment,"auxiliary commitment","pivot","terminal principal equation"},
                             {Fault::residual,"invalid fresh residual","field","epoch sumcheck recurrence"},
                             {Fault::pending,"fabricated pending increment","pivot","terminal principal equation"},
                             {Fault::continuity,"broken state link","collection","state continuity"},
                             {Fault::response,"altered pivot response","pivot","pivot binding equation"}}) {
      std::cout << "CASE " << test.name << '\n';
      try {
        run_case(master,test.fault);
      } catch (const VerifierRejection& e) {
        if (std::string_view(e.stage)!=test.stage || std::string_view(e.what())!=test.rejection) throw;
        std::cout << "REJECT " << test.name << " (" << e.stage << ": " << e.what() << ")\n";
        continue;
      }
      sonata::check(!*test.rejection,"mutation unexpectedly accepted");
      std::cout << "ACCEPT " << test.name << '\n';
    }
    std::cout << "PASS: small live verifier accepted the honest history and rejected all mutations\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
}
