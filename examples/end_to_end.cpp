#include "challenge_fixture.hpp"
#include "five_layer.hpp"
#include "window_input.hpp"
#include <cstdlib>
#include <exception>
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
  if (argc != 2 || std::string_view(argv[1]) != "--run-full" ||
      !pipe_mode || std::string_view(pipe_mode) != "1") {
    std::cerr << "Usage: python3 examples/run_e2e.py ./build/sonata_e2e --run-full\n";
    return 2;
  }
  try {
    const auto crs = example_crs();
    std::mt19937_64 field_rng(0x534f4e41), backend_rng(0x54414532);
    std::cerr << "Preparing the eight-transition checkpoint (994 fresh instances)\n";
    auto input = sonata::example::make_window_input(crs, field_rng);
    sonata::example::run_five_layers(
        std::move(input.system), std::move(input.witness), crs,
        [&backend_rng] { return backend_rng(); }, sonata::example::draw_challenges);
    std::cout << "PASS: checkpoint verified through five layers and altered pivot rejected\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
