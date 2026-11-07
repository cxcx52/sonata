#include "sonata/epoch_envelope.hpp"
#include <iostream>

using namespace sonata;

int main() {
  // Active coefficients 0..65, ignored 66..128, mandatory zero 129..191.
  // Both cuts fall inside a ring.
  Config cfg{1,3,0,129,192,0,0,0,0};
  System active;
  active.raw = 66;
  active.y = {17, 1};
  active.apply = [](const RV &v) -> SV {
    check(v.size() == 2, "active ring count");
    return {add(cv(v,1), mul(3,cv(v,64))), cv(v,65)};
  };
  active.adj = [](const SV &a) {
    check(a.size() == 2, "active rows");
    RV z(2);
    put(z,1,a[0]); put(z,64,mul(3,a[0])); put(z,65,a[1]);
    for (auto &ring : z) ring = dual(ring);
    return z;
  };
  System root = align(epoch::embed_prefix(active,cfg,70),cfg);
  RV w(3);
  put(w,1,2); put(w,64,5); put(w,65,1);
  check(root.apply(w) == root.y, "honest zero extension");
  put(w,66,4); put(w,128,7);
  check(root.apply(w) == root.y, "different ignored tails, same active label");
  put(w,129,9);
  check(root.apply(w) != root.y && root.apply(w).at(7)==9, "original align padding retained");
  for (size_t i=0; i<192; ++i) put(w,i,(i*i+19)%q);
  auto out=root.apply(w);
  for (size_t i=0; i<70; ++i) {
    SV a(70); a[i]=1;
    check(ctDot(root.adj(a),w)==out[i], "padded row transpose");
  }
  std::cout << "PASS: ignored-tail projection, original align and all 70 adjoint rows\n";
}
