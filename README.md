![Sonata](assets/sonata.png)

# Sonata

Sonata folds an accumulator and a batch of fresh lattice instances into a new accumulator. Each transition checks the old and new state equations and records a weighted sum of the fresh residuals in a public field element. In a sequence of authenticated transitions, a final zero-pending decision checks the outstanding residuals.

The repository provides C++17 prover and verifier code, a Python short-challenge sampler and audited five-layer profiles for batches of 16, 64, 256 and 1,024 fresh instances. The knowledge analysis uses labelled two-point extraction against a classically rewindable prover.

## Five-layer results

Each ordinary transition starts from committed inputs, draws fresh verifier challenges, and completes four recursive reductions and the terminal pivot. The verifier reconstructs all five layers from public inputs and the decoded messages. Times are medians of three runs per arity.

**Measurement environment.** The shared server used AMD EPYC 7742 processors (x86-64), Ubuntu 24.04.4 LTS and Linux `5.15.0-161-generic`. The Release build used GCC 13.3.0, CMake 3.28.3 and `-O3 -DNDEBUG`. The verifier driver used Python 3.12.3. The container had a 200 GiB memory limit and a CPU time quota equivalent to 40 logical CPUs. Within each arity, the three trials ran concurrently, one each on logical CPUs 4, 36 and 68, on distinct physical cores with separate L3 caches.

| Fresh inputs | Prover (s) | Verifier (s) | Wall (s) | Bidirectional transcript (B) | Peak RSS (MiB) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 19.737 | 18.071 | 37.912 | 143,165 | 86.65 |
| 64 | 34.227 | 31.068 | 65.544 | 149,246 | 167.50 |
| 256 | 83.797 | 74.946 | 159.603 | 159,507 | 523.89 |
| 1,024 | 258.324 | 226.150 | 487.709 | 170,361 | 2,014.43 |

Bytes count the two encoded message streams with input commitments and CRS already shared. Peak RSS is the largest single-process high-water mark across the three trials. At 1,024 inputs the actual root commitment matrix has 16 rows and 4,733 columns, and the shared packed CRS occupies 29,079,552 B. [Full tables](benchmarks/results/experiment-table.md) · [Individual trials](benchmarks/results/trials.csv) · [Reproduction commands](#four-arity-ordinary-transition-experiment)

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

**Measurement environment.** The component benchmarks ran on an Intel Core i9-13980HX under Ubuntu 24.04 LTS on WSL2, kernel `6.18.33.2-microsoft-standard-WSL2`. The C++ kernels used GCC 13.3.0 with `-O3 -DNDEBUG -std=gnu++17`, and the short-challenge sampler used CPython 3.12.3. Each module ran sequentially in five fresh processes pinned to logical CPU 4, without a warm-up. Time is the median for the listed batch. Memory is the maximum whole-process peak RSS, including input preparation and runtime overhead.

| Operation | Measured input | Time (ms) | Peak RSS (MiB) |
| --- | --- | ---: | ---: |
| Canonical encoding into the opening | 66,064 field scalars | 4.073 | 31.28 |
| Square Sumcheck prover on a prepared table | 44,386 entries, 16 rounds | 4.502 | 4.53 |
| Seed expansion and projection of one vector | 384 × 44,160 | 799.013 | 19.53 |
| Short-challenge sampling with exact spectral filtering | 21 rings of degree 64, prepared rank table | 11.704 | 19.31 |

Build with `cmake --build build --target sonata_bench` and run `python3 benchmarks/run.py ./build/sonata_bench` on Linux with GNU time. The [benchmark driver](benchmarks/run.py) emits individual samples and the machine configuration as JSON.

### Manual checkpoint example

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

The window run requires exit status zero and the final line `PASS: checkpoint verified through five layers and altered pivot rejected`. To record the wall time and peak RSS for this fixture, run `python3 benchmarks/full_profile.py window ./build/sonata_e2e window-results` with a new output directory.

### Four-arity ordinary-transition experiment

The four profiles share a 128-bit classical attack-cost target. All six generic Core-SVP estimates must reach 144 under classical pricing and 136 under quantum pricing. The audit also checks single-CRT-component attacks and requires a statistical error of at most $2^{-130}$ over 16,384 issued attempts, including abandoned attempts. Quantum pricing describes the cost of lattice reduction. The classical knowledge theorem retains a global extraction cutoff and six time-bounded Module-SIS advantages. [Parameter ledger](parameters/sweep-profiles.json) · [Attack estimates and MATZOV fit ranges](benchmarks/results/attacks.csv)

The batch driver uses fixed public input fixtures with one nonzero fresh residual and tests rejection of an altered pivot commitment in every run. Use Python 3.11 or later for the parameter audit and benchmark scripts. For the attack estimates, use SageMath 10.9 and a clean Lattice Estimator checkout at the pinned revision. Replace the CPU IDs below with three available cores with distinct last-level caches and choose fresh output directories. [Batch driver](benchmarks/run_four_arities.py) · [Pinned Lattice Estimator](https://github.com/malb/lattice-estimator/tree/53da5982597709ba0fdf94ea37a84d822310fd84)

```sh
python3 parameters/sweep_profiles.py --check
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target sonata_sweep
python3 -B benchmarks/run_four_arities.py build/sonata_sweep runs --cpus 4 36 68
sage -python parameters/estimate_sis.py /path/to/lattice-estimator parameters/sweep-profiles.json > sis-sweep.json
python3 -B benchmarks/summarize_four_arities.py runs sis-sweep.json tables
```

The [summary script](benchmarks/summarize_four_arities.py) independently decodes each stream, checks lengths against the issued retry sequence, and writes per-trial and aggregate results, packed CRS sizes, and 24 attack-cost rows with MATZOV fit-range flags. Prover and verifier times exclude shared CRS and input preparation. The verifier column includes issuing the coins and checking the decoded transcript without the witness. The table also reports public replay alone and the complete wall time. Wire lengths include padding and retries, and peak RSS is the largest individual-process high-water mark across the three trials.

## Code

| Component | Entry point |
| --- | --- |
| Field and negacyclic ring arithmetic | [field.hpp](include/sonata/field.hpp) · [ring.hpp](include/sonata/ring.hpp) |
| Pending transition and Sumcheck endpoint | [frontend.hpp](include/sonata/frontend.hpp) |
| Seeded projection, four-coordinate aggregation, child and pivot checks | [backend.hpp](include/sonata/backend.hpp) |
| Short-challenge sampler and codec | [ball_challenges.py](extraction/ball_challenges.py) |

## Origins

Building on CauchyFold and an unpublished internal study of knowledge extraction, Sonata combines pending-residual folding with a two-point extractor.

The original CauchyFold code, parameters, paper citation and MIT license are preserved together. To reproduce its results, start in [archive/cauchyfold](archive/cauchyfold/README.md).

I would be grateful to anyone who helps carry out the remaining experiments or shares feedback on Sonata.
