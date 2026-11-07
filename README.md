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

## Checkpoint batching

A verifier can freeze a window of consecutive transitions and submit one root statement to the backend. The [field session](include/sonata/epoch_session.hpp) takes its own copy of the window before issuing the batching challenges. The outer verifier supplies fresh randomness, shares the attempt budget across retries, and accepts the checkpoint only after the backend accepts the statement.

For eight transitions with **994 fresh instances in total**, the window uses 22 Sumcheck rounds. Certifying the same eight transitions separately uses 152 rounds. The [capacity calculator](parameters/epoch.py) checks which window sizes fit the existing five-layer backend. Run `python parameters/epoch.py --check` to reproduce the [capacity ledger](parameters/epoch-capacity.json).

## Build and test

The C++ kernels build on Linux or WSL with CMake 3.16+ and GCC or Clang with native 128-bit integer support.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake -E chdir build ctest --output-on-failure
```

The [four-input example](tests/cpp/pending_child.cpp) follows one pending transition through its four-coordinate child relation to the terminal pivot. Its fixtures use deterministic verifier challenges. Applications supply the common reference string and the verifier's challenges in protocol order.

To reproduce the eight-transition comparison, run `./build/sonata_epoch --target` for the window and `./build/sonata_epoch --target-sequential` for the eight separate certifications. Both commands check the frontend through its root equations. Use `--target-direct` to exercise the window variant that eliminates the pending equations when consecutive public pending values agree.

## Code

| Component | Entry point |
| --- | --- |
| Field and negacyclic ring arithmetic | [field.hpp](include/sonata/field.hpp) · [ring.hpp](include/sonata/ring.hpp) |
| Pending transition and Sumcheck endpoint | [frontend.hpp](include/sonata/frontend.hpp) |
| Seeded projection, four-coordinate aggregation, child and pivot checks | [backend.hpp](include/sonata/backend.hpp) |
| Short-challenge sampler and codec | [ball_challenges.py](extraction/ball_challenges.py) |
| Local slope recovery and short-kernel comparison from three authenticated records | [coordinate.py](extraction/coordinate.py) |

Run `python -m unittest discover -s tests/python` to check the Python coordinate routines.

## Origins

Building on CauchyFold and an unpublished internal study of knowledge extraction, Sonata combines pending-residual folding with a two-point extractor.

The original CauchyFold code, parameters, paper citation and MIT license are preserved together. To reproduce its results, start in [archive/cauchyfold](archive/cauchyfold/README.md).

I would be grateful to anyone who helps carry out the remaining experiments or shares feedback on Sonata.
