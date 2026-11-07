![Sonata](assets/sonata.png)

# Sonata

Sonata folds an accumulator and a batch of fresh lattice instances into a new accumulator. Each transition checks the old and new state equations and records a weighted sum of the fresh residuals in a public field element. In a sequence of authenticated transitions, a final zero-pending decision checks the outstanding residuals.

The repository provides C++17 kernels for the interactive relations, Python routines for short challenges and local extraction, and the concrete parameters for a batch of 1,024 fresh instances.

## Reproduce the parameters

From the repository root, run the following with Python 3.11 or later. The calculator uses only the standard library.

```sh
python parameters/compute.py --points 2 --check
python parameters/compute.py --points 3 --check
```

Each command checks the full five-layer calculation against its committed JSON file. Omit `--check` to print the dimensions, norm bounds and exact statistical error terms.

| Five-layer profile | Two-point extraction | Three-point extraction |
| --- | ---: | ---: |
| Zero-pending statistical knowledge error, −log₂(bound) | 144.18 | 144.18 |
| Upper-bound coefficient on expected terminal calls | 2,276,736 | 53,717,664 |
| Parameters | [Two-point ledger](parameters/ledger-two-point.json) | [Three-point ledger](parameters/ledger.json) |

Both variants use the same interactive messages. The terminal-call coefficient bounds work in the extractor. Concrete computational security also depends on the Module-SIS instances and the global extraction time limit.

## Build and test

The C++ kernels build on Linux or WSL with CMake 3.16+ and GCC or Clang with native 128-bit integer support.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake -E chdir build ctest --output-on-failure
```

The [four-input example](tests/cpp/pending_child.cpp) follows one pending transition through its four-coordinate child relation to the terminal pivot. Its fixtures use deterministic verifier challenges. Applications supply the common reference string and the verifier's challenges in protocol order.

## Code

| Component | Entry point |
| --- | --- |
| Field and negacyclic ring arithmetic | [field.hpp](include/sonata/field.hpp) · [ring.hpp](include/sonata/ring.hpp) |
| Pending transition and Sumcheck endpoint | [frontend.hpp](include/sonata/frontend.hpp) |
| Seeded projection, four-coordinate aggregation, child and pivot checks | [backend.hpp](include/sonata/backend.hpp) |
| Short-challenge sampler and codec | [ball_challenges.py](extraction/ball_challenges.py) |
| Local slope recovery and short-kernel comparison from three authenticated records | [coordinate.py](extraction/coordinate.py) |

Run `python -m unittest discover -s tests/python` to check the Python coordinate routines.

## CauchyFold

Sonata builds on CauchyFold and Palimpsest. The original CauchyFold code, parameters, paper citation and MIT license are preserved together. To reproduce its results, start in [archive/cauchyfold](archive/cauchyfold/README.md).
