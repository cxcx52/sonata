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

Run `python parameters/audit_profile.py` for an independent reconstruction of both ledgers, including the short-challenge support, balanced-digit norm bounds, six representative Module-SIS instances, and per-layer communication. The packed backend messages total 166,136 bytes with no retries. Settling the eight-step window below adds 3,192 bytes for the field exchange and auxiliary commitment when the verifier already holds the public log. These are encoding budgets.

The [attack estimator driver](parameters/estimate_sis.py) checks all six representatives against the nonzero-modulo-$q$ kernel condition and their Euclidean radii. It runs the coefficient-lattice model of the pinned [Lattice Estimator](https://github.com/malb/lattice-estimator/tree/53da5982597709ba0fdf94ea37a84d822310fd84) and reports the BKZ block sizes, cost-model sensitivity, and a separate single-CRT-component screen. Use SageMath 10.9 and a clean checkout at that commit to reproduce the estimates.

```sh
sage -python parameters/estimate_sis.py /path/to/lattice-estimator parameters/ledger-two-point.json
```

The smallest classical Core-SVP estimate among the six generic models is 142.8 in log₂ work units ($A_1$, block size 489). Bounding the protocol failure probability also requires six time-bounded success probabilities and the history extractor's cutoff loss.

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

### Small live verifier

The [live example](examples/interactive.cpp) certifies two consecutive steps containing one and two fresh instances. It runs the 16-round field Sumcheck, one recursive reduction and the terminal pivot. The verifier fixes the incoming commitments before drawing the fold challenges, generates each later challenge after its preceding message, and uses fresh operating-system randomness. The prover and verifier run in one process, with a Python worker supplying filtered short challenges to the verifier.

```sh
cmake --build build --target sonata_interactive
python3 examples/run_interactive.py ./build/sonata_interactive
```

The run checks acceptance of an honest window and rejection of an altered auxiliary commitment, a committed source with a false residual, a fabricated intermediate pending value, a broken state link and an altered terminal response.

### Component benchmarks

Each value below comes from five fresh processes on one pinned logical CPU of an Intel Core i9-13980HX under Ubuntu 24.04 / WSL2. Time is the median for the listed batch. Memory is the maximum whole-process peak RSS, including input preparation and runtime overhead.

| Operation | Measured input | Time (ms) | Peak RSS (MiB) |
| --- | --- | ---: | ---: |
| Canonical encoding into the opening | 66,064 field scalars | 4.073 | 31.28 |
| Square Sumcheck prover on a prepared table | 44,386 entries, 16 rounds | 4.502 | 4.53 |
| Seed expansion and projection of one vector | 384 × 44,160 | 799.013 | 19.53 |
| Short-challenge sampling with exact spectral filtering | 21 rings of degree 64, prepared rank table | 11.704 | 19.31 |

Build with `cmake --build build --target sonata_bench` and run `python3 benchmarks/run.py ./build/sonata_bench` on Linux with GNU time. The [benchmark driver](benchmarks/run.py) emits individual samples and the machine configuration as JSON.

### Manual end-to-end example

The [full-profile example](examples/end_to_end.cpp) starts from a trusted fixture state and constructs eight linked transitions containing 994 fresh instances. The initial and final pending values are zero. It freezes the window, runs the field Sumcheck, and passes the resulting root and opening through four recursive reductions and the terminal pivot. A separate verifier reconstructs the backend relations from the public root, messages, and its saved challenges. It checks the honest pivot, then changes one commitment coefficient and checks for rejection. The fixture uses deterministic verifier coins for reproducibility.

Compile the example on 64-bit Linux or WSL with the following commands.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target sonata_e2e
```

The target is excluded from the default build and CTest. The root projection alone stores 384 × 3,635,520 one-byte coefficients (1.30 GiB). Allow additional memory for the witness, commitments, and aggregation. The output includes the public replay of all five layers and one deliberately altered terminal commitment.

On a machine reserved for the full run, start the [Python driver](examples/run_e2e.py) from the repository root. It supplies the spectrally filtered ring challenges to the C++ executable.

```sh
python3 examples/run_e2e.py ./build/sonata_e2e --run-full
```

For one ordinary transition with 1,024 fresh instances, add `--transition`. This mode draws fresh verifier coins and uses a source with one nonzero residual to exercise the pending update. To collect wall time, the largest child process RSS high-water mark, and process-tree RSS samples on Linux, use a new output directory for each run.

```sh
python3 benchmarks/full_profile.py transition ./build/sonata_e2e transition-results
python3 benchmarks/full_profile.py window ./build/sonata_e2e window-results
```

On an AMD EPYC 7742 with a Release build, one complete ordinary run took 523.79 s (2.062 GiB largest-process peak RSS); the deterministic checkpoint run took 524.78 s (1.632 GiB). These are separate local-message harness measurements covering proof creation, public verification and a negative pivot check.

The window run requires exit status zero and the final line `PASS: checkpoint verified through five layers and altered pivot rejected`. This exercises one honest interactive execution and one tampered terminal message. Knowledge extraction has a separate parameter ledger and local tests.

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
