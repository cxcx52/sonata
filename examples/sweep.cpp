#include "challenge_fixture.hpp"
#include "sweep_profiles.hpp"
#include "sweep_wire.hpp"
#include "transition_input.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace {
using namespace sonata;
using namespace sonata::example;
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Proof = std::array<LayerMessages, 5>;
using Tape = std::array<LayerCoins, 5>;

double since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

SweepProfile select_profile(size_t k, std::string_view mode) {
  if (mode == "reference") {
    check(k == 1024, "reference profile requires k=1024");
    return {1024, 44352, layers};
  }
  check(mode == "optimized", "unknown profile family");
  for (const auto& p : sweep_profiles)
    if (p.k == k) return p;
  throw std::runtime_error("unsupported arity");
}

void save(const fs::path& file, const Bytes& bytes) {
  std::ofstream out(file, std::ios::binary);
  out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
  out.close();
  check(bool(out), "could not write binary transcript");
}

Bytes read(const fs::path& file) {
  std::ifstream in(file, std::ios::binary);
  check(bool(in), "could not open transcript");
  Bytes bytes(fs::file_size(file));
  in.read(reinterpret_cast<char*>(bytes.data()), bytes.size());
  check(bool(in) && in.peek() == EOF, "incomplete transcript file");
  return bytes;
}

void attempt_flags(WireWriter& writer, size_t retries) {
  for (size_t j = 0; j < retries; ++j) writer.flag(false);
  writer.flag(true);
}

WireWriter encode_prover(const PublicFront& f, const Proof& proof, const Tape& tape,
                         const LayerPath& path) {
  WireWriter writer;
  for (const E& v : f.pub.back()) writer.field(v);
  writer.rings(f.cs.back());
  writer.field(f.pending);
  writer.rings(f.caux);
  for (size_t j = 0; j < f.field.ch.size(); ++j) {
    writer.field(f.field.g0[j]);
    writer.field(f.field.g2[j]);
  }
  writer.field(f.field.terminal);
  for (size_t i = 0; i < path.size(); ++i) {
    const auto& cfg = path[i];
    const auto& msg = proof[i];
    if (cfg.rho) writer.rings(msg.u1);
    else {
      writer.rings(msg.pivot.pivot);
      writer.rings(msg.pivot.tothers);
    }
    attempt_flags(writer, tape[i].skipped_seeds.size());
    writer.bounded(msg.image, projrows*cfg.S);
    if (cfg.rho) writer.rings(msg.u2);
    else writer.triangle(msg.pivot.h, cfg.s);
    attempt_flags(writer, tape[i].skipped_challenges.size());
    if (!cfg.rho) {
      SV coefficients;
      append(coefficients, msg.pivot.z);
      writer.bounded(coefficients, cfg.G);
    }
  }
  return writer;
}

struct FieldCoins {
  E alpha, gamma, rho, eta;
  std::vector<E> rounds;
};

WireWriter encode_verifier(const FieldCoins& f, const Tape& tape) {
  WireWriter writer;
  for (const auto& v : {f.alpha, f.gamma, f.rho, f.eta}) writer.field(v);
  for (const auto& v : f.rounds) writer.field(v);
  for (const auto& coins : tape) {
    for (const auto& v : coins.skipped_seeds) writer.seed(v);
    writer.seed(coins.projection_seed);
    writer.field(coins.zeta);
    for (const auto& v : coins.skipped_challenges) writer.challenges(v);
    writer.challenges(coins.challenge);
  }
  return writer;
}

// The verifier starts with the original public inputs, its retained coins and
// the two encoded streams. All output-side records are reconstructed below.
void replay(PublicFront f, const FieldCoins& fields, const Tape& issued,
            const std::shared_ptr<const Master>& crs, const SweepProfile& profile,
            const Bytes& p_bytes, const Bytes& v_bytes, double& tamper_time) {
  WireReader p(p_bytes), v(v_bytes);
  auto receive_coin = [&](const E& original) {
    E decoded = v.field();
    check(decoded == original, "field challenge differs from the issued value");
    return decoded;
  };
  f.alpha = receive_coin(fields.alpha);
  f.gamma = receive_coin(fields.gamma);
  std::array<E, 5> output;
  for (auto& e : output) e = p.field();
  f.pub.push_back(output);
  f.cs.push_back(p.rings(4));
  f.pending = p.field();
  f.caux = p.rings(4);
  f.rho = receive_coin(fields.rho);
  f.eta = receive_coin(fields.eta);
  for (const E& r : fields.rounds) {
    f.field.g0.push_back(p.field());
    f.field.g2.push_back(p.field());
    f.field.ch.push_back(receive_coin(r));
  }
  f.field.terminal = p.field();
  auto root = align(initialsystem(f, crs), profile.layers.front());
  Proof proof{};
  Tape accepted{};
  for (size_t i = 0; i < profile.layers.size(); ++i) {
    const auto& cfg = profile.layers[i];
    auto& msg = proof[i];
    const auto& coins = issued[i];
    if (cfg.rho) msg.u1 = p.rings(6);
    else {
      msg.pivot.pivot = p.rings(2);
      msg.pivot.tothers = p.rings(16*(cfg.s-1));
    }
    const size_t seed_bytes = (cfg.b + cfg.h*(3*cfg.b-1) + 7)/8;
    check(coins.skipped_seeds.size() < attempt_limit, "excess projection attempts");
    for (const auto& seed : coins.skipped_seeds)
      check(v.seed(seed_bytes) == seed && !p.flag(), "projection retry transcript");
    accepted[i].projection_seed = v.seed(seed_bytes);
    check(accepted[i].projection_seed == coins.projection_seed && p.flag(), "projection submission");
    msg.image = p.bounded(projrows, projrows*cfg.S);
    accepted[i].zeta = receive_coin(coins.zeta);
    if (cfg.rho) msg.u2 = p.rings(6);
    else msg.pivot.h = p.triangle(cfg.s);
    check(coins.skipped_challenges.size() < attempt_limit, "excess response attempts");
    for (const auto& c : coins.skipped_challenges)
      check(v.challenges(cfg.s) == c && !p.flag(), "short challenge retry transcript");
    accepted[i].challenge = v.challenges(cfg.s);
    check(accepted[i].challenge == coins.challenge && p.flag(), "short challenge submission");
    if (!cfg.rho) msg.pivot.z = chunks(p.bounded(degree*cfg.n, cfg.G), 0, cfg.n);
  }
  check(p.done() && v.done(), "unexpected trailing transcript data");
  verify_layers(std::move(root), crs, proof, accepted, profile.layers, &tamper_time);
}

void run(const SweepProfile& profile, const fs::path& dir, bool reference) {
  check(fs::is_directory(dir), "driver must supply a fresh result directory");
  const auto start = Clock::now();
  auto crs = std::make_shared<Master>();
  std::mt19937_64 crs_rng(79);
  crs->data.resize(profile.crs_cells);
  for (auto& ring : crs->data) ring = fixture::randR(crs_rng);
  std::ifstream entropy("/dev/urandom", std::ios::binary);
  auto word = [&] {
    U value;
    check(bool(entropy.read(reinterpret_cast<char*>(&value), sizeof(value))), "verifier entropy unavailable");
    return value;
  };

  PublicFront common{};
  FieldCoins field_tape;
  Tape tape{};
  WireWriter prover_wire;
  double setup = 0, coins = 0, proving = 0;
  {
    auto input = prepare_ordinary_inputs(static_cast<int>(profile.k), crs);
    common = input.pub;
    setup = since(start);
    const auto begin = Clock::now();
    auto [root, w] = finish_ordinary_transition(input, crs, word, profile.layers[0], &coins);
    field_tape = {input.pub.alpha, input.pub.gamma, input.pub.rho, input.pub.eta, input.pub.field.ch};
    Proof proof{};
    prove_layers(std::move(root), std::move(w), crs, word, draw_challenges, proof, tape,
                 profile.layers, &coins);
    prover_wire = encode_prover(input.pub, proof, tape, profile.layers);
    proving = since(begin) - coins;
  }
  const auto encoding_start = Clock::now();
  auto verifier_wire = encode_verifier(field_tape, tape);
  const double verifier_encode = since(encoding_start);
  const fs::path p_file = dir/"p-to-v.bin", v_file = dir/"v-to-p.bin";
  save(p_file, prover_wire.data);
  save(v_file, verifier_wire.data);
  const auto p_bytes = read(p_file), v_bytes = read(v_file);
  const auto verify_start = Clock::now();
  double tamper = 0;
  replay(std::move(common), field_tape, tape, crs, profile, p_bytes, v_bytes, tamper);
  const double verification = since(verify_start) - tamper;
  std::ofstream report(dir/"metrics.json");
  report << std::setprecision(12)
         << "{\n  \"k\": " << profile.k << ",\n  \"reference\": " << (reference ? "true":"false")
         << ",\n  \"prover_seconds\": " << proving
         << ",\n  \"verifier_seconds\": " << verification + coins + verifier_encode
         << ",\n  \"verifier_coin_seconds\": " << coins
         << ",\n  \"verifier_encode_seconds\": " << verifier_encode
         << ",\n  \"verify_replay_seconds\": " << verification
         << ",\n  \"setup_seconds\": " << setup
         << ",\n  \"tamper_seconds\": " << tamper
         << ",\n  \"prover_bytes\": " << p_bytes.size()
         << ",\n  \"verifier_bytes\": " << v_bytes.size();
  for (bool seeds : {true, false}) {
    report << ",\n  \"" << (seeds ? "projection_attempts" : "short_attempts") << "\": [";
    for (size_t i = 0; i < tape.size(); ++i) {
      if (i) report << ", ";
      report << 1 + (seeds ? tape[i].skipped_seeds.size() : tape[i].skipped_challenges.size());
    }
    report << ']';
  }
  report << "\n}\n";
  report.close();
  check(bool(report), "could not write measurements");
}
} // namespace

int main(int argc, char** argv) {
  try {
    const char* pipe = std::getenv("SONATA_E2E_CHALLENGE_PIPE");
    sonata::check(argc == 4 && pipe && std::string_view(pipe) == "1", "use examples/run_sweep.py");
    const auto profile = select_profile(std::stoul(argv[1]), argv[3]);
    run(profile, argv[2], std::string_view(argv[3]) == "reference");
    std::cout << "PASS: verified five layers from the encoded transition and rejected the altered pivot\n";
  } catch (const std::exception& e) {
    std::cerr << "FAIL: " << e.what() << '\n';
    return 1;
  }
}
