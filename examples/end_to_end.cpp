#include "challenge_fixture.hpp"
#include "five_layer.hpp"
#include "transition_input.hpp"
#include "window_input.hpp"
#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <random>
#include <string_view>

namespace {

std::shared_ptr<const sonata::Master> example_crs() {
  std::mt19937_64 rng(79);
  auto crs = std::make_shared<sonata::Master>();
  crs->data.resize(44352);
  for (auto& cell : crs->data)
    cell = sonata::fixture::randR(rng);
  return crs;
}

} // namespace

int main(int argc, char** argv) {
  const char* pipe_mode = std::getenv("SONATA_E2E_CHALLENGE_PIPE");
  const std::string_view mode = argc == 2 ? argv[1] : "";
  const bool transition = mode == "--transition";
  if ((!transition && mode != "--run-full") ||
      !pipe_mode || std::string_view(pipe_mode) != "1") {
    std::cerr << "Usage: python3 examples/run_e2e.py ./build/sonata_e2e --run-full [--transition]\n";
    return 2;
  }
  try {
    const auto start = std::chrono::steady_clock::now();
    const auto crs = example_crs();
    const auto crs_ready = std::chrono::steady_clock::now();
    std::mt19937_64 field_rng(0x534f4e41), backend_rng(0x54414532);
    std::ifstream entropy;
    sonata::example::WordSource word = [&] { return backend_rng(); };
    sonata::example::WindowInput input;
    if (transition) {
      std::cerr << "Preparing one ordinary transition (1024 fresh instances)\n";
      entropy.open("/dev/urandom", std::ios::binary);
      word = [&] {
        sonata::U value;
        if (!entropy.read(reinterpret_cast<char*>(&value), sizeof(value)))
          throw std::runtime_error("verifier entropy unavailable");
        return value;
      };
      auto [root, witness] = sonata::example::ordinary_transition(crs, word);
      input = {std::move(root), std::move(witness)};
    } else {
      std::cerr << "Preparing the eight-transition checkpoint (994 fresh instances)\n";
      input = sonata::example::make_window_input(crs, field_rng);
    }
    const auto front_ready = std::chrono::steady_clock::now();
    std::cerr << "crs_seconds=" << std::chrono::duration<double>(crs_ready-start).count()
              << " frontend_seconds=" << std::chrono::duration<double>(front_ready-crs_ready).count() << '\n';
    sonata::example::run_five_layers(
        std::move(input.system), std::move(input.witness), crs,
        word, sonata::example::draw_challenges);
    std::cout << (transition ? "PASS: ordinary transition" : "PASS: checkpoint")
              << " verified through five layers and altered pivot rejected\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
