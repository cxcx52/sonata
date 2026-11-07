#pragma once
#include "sonata/backend.hpp"
#include <random>

// Single terminal layer to check the small epoch's real backend interface.
// This fixture has different dimensions from the target five-layer profile.
inline void check_epoch_pivot(const sonata::System &active,sonata::RV w,
                              const std::shared_ptr<sonata::Master>& master) {
  using namespace sonata;
  const size_t s=16, n=(active.raw+64*s-1)/(64*s);
  const U S=47*((active.raw-1)/55)+1;
  Config cfg{s,n,0,active.raw,64*s*n,S,256*S,1174,16};
  auto aligned=align(active,cfg);
  w.resize(s*n);
  Key A{master,16,n}, pivot{master,2,256};
  TerminalResponse proof;
  for(size_t j=0;j<s;++j){
    auto image=A.apply(slice(w,j*n,n));
    if(!j) proof.pivot=pivot.apply(digits(image,8,16));
    else proof.tothers.insert(proof.tothers.end(),image.begin(),image.end());
  }
  size_t bits=cfg.b+cfg.h*(3*cfg.b-1);
  Bytes seed((bits+7)/8);
  std::mt19937_64 rng(1941);
  for(auto& b:seed) b=static_cast<uint8_t>(rng());
  if(bits%8) seed.back()&=(1u<<(bits%8))-1;
  auto projection=expand_projection(cfg,seed);
  auto p=project(projection,w);
  check(proj_norm(p,S),"small epoch projection");
  auto ag=aggregate(aligned,projection,p,E{11,8,0,2},one);
  proof.h=symmetric(ag,w,s,n);
  for(size_t e=0;e<4;++e) proof.h[e*(s*(s+1)/2)][0]=0;
  RV c(s);
  for(size_t j=0;j<s;++j) c[j][j]=j%2?q-1:1;
  proof.z=response(w,c,n);
  verifyterminal(cfg,A,pivot,ag,proof,c);
  proof.h[tri(s,0,1)][0]=add(proof.h[tri(s,0,1)][0],1);
  bool rejected=false;
  try { verifyterminal(cfg,A,pivot,ag,proof,c); }
  catch(const std::runtime_error&) { rejected=true; }
  check(rejected,"small epoch pivot accepted altered principal table");
  std::cout<<"PASS: two-step epoch opening through seeded four-lane pivot and mutation\n";
}
