#include "../../examples/sweep_wire.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace sonata;
using namespace sonata::example;

template <class F> static void rejects(F decode) {
  try { decode(); } catch (const std::runtime_error &) { return; }
  throw std::logic_error("malformed transcript accepted");
}

int main() {
  WireWriter p;
  SV signed_values{0,q-1,1,q-2,3};
  p.bounded(signed_values,30);
  check(p.data==Bytes({0x59,0xaf,0}),"Rice fixture");
  check(wire_detail::rice_layout(384,384ULL*3105761).bytes==656,"projection budget");
  check(wire_detail::rice_layout(29696,113273922816ULL).bytes==51520,"response budget");
  Ring r{}; r[0]=q-1; r[1]=2; r[2]=q-3;
  WireWriter v;
  v.challenges(RV{r});
  std::ostringstream hex;
  for (auto byte:v.data) hex << std::hex << std::setw(2) << std::setfill('0') << unsigned(byte);
  check(hex.str()=="c3238173fe05056edd58c2b7ed77329d96bdfe35","Python BallIndex rank");
  const E field{0,1,255,q-1};
  RV triangle(12, r);
  for (size_t j=0;j<12;j+=3) triangle[j][0]=0;
  v.field(field); v.rings(RV{r}); v.triangle(triangle,2);
  v.seed(Bytes{0x80,0x01}); v.flag(false); v.flag(true);
  WireReader vr(v.data), pr(p.data);
  check(vr.challenges(1)==RV{r} && vr.field()==field && vr.rings(1)==RV{r} &&
        vr.triangle(2)==triangle && vr.seed(2)==Bytes({0x80,0x01}) &&
        !vr.flag() && vr.flag() && vr.done(),"typed wire roundtrip");
  check(pr.bounded(5,30)==signed_values && pr.done(),"Rice roundtrip");
  auto bad=p.data; bad.back()=0x80;
  rejects([&]{ WireReader(bad).bounded(5,30); });
  bad.pop_back();
  rejects([&]{ WireReader(bad).bounded(5,30); });
  Bytes invalid(24,0xff);
  rejects([&]{ WireReader(invalid).field(); });
  rejects([&]{ WireReader(invalid).challenges(1); });
  rejects([&]{ WireReader(invalid).flag(); });
  std::cout<<"PASS packed messages, independent rank and Rice fixtures, malformed input rejection\n";
}
