# Collie

Collie is a confidence-level limit and profile-likelihood toolkit for binned
particle-physics analyses. It separates ROOT file I/O from statistical model
evaluation and supports deterministic, parallel pseudo-experiments.

## Requirements

- A C++17 compiler
- CMake 3.20 or newer
- ROOT 6.20 or newer, including Minuit2

Collie no longer builds private copies of CLHEP or Fortran MINUIT. Random
generation uses the C++ standard library, and profile fits use ROOT Minuit2.

## Build and test

Load a ROOT environment so `root-config` is available, then run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Install to a local prefix with:

```sh
cmake --install build --prefix "$PWD/install"
```

The original Make workflow remains available:

```sh
make -j
source setup_Collie.sh
```

## Formatting

The repository uses a checked-in `clang-format` policy for C++ and ROOT macros.
Format the source tree or verify it without changing files with:

```sh
./scripts/format.sh
./scripts/format.sh --check
```

## Unified command line

The `collie` executable replaces the collection of one-off drivers for new
workflows:

```sh
collie inspect input.root
collie validate --strict input.root
collie limit --method fast --seed 1234 --workers 8 \
  --output limits.root input.root
collie fit --point 300 --output fit.root input.root
collie combine --output combined.root first.root second.root
```

Run `collie COMMAND --help` for all options. Outputs are written atomically, so a
failed run does not leave a partial ROOT file. Each produced file contains a
`collie_run_manifest` with the master seed, worker count, algorithm settings,
ROOT version, and input checksums.

## Reproducibility and parallel execution

Stochastic work is controlled by `collie::RunContext`. Set it through CLI
options or environment variables:

```sh
export COLLIE_RANDOM_SEED=1234
export COLLIE_NUM_THREADS=8
```

Pseudo-experiment seeds are derived from the master seed and trial index.
Consequently, the parallel `CLfast` result is identical for one worker or many
workers; scheduling does not change the generated experiment stream. Exact
bitwise replay assumes the same standard-library implementation because the
C++ standard does not prescribe identical distribution algorithms across
runtimes.

## Library structure

- `CollieIO` reads and writes the established ROOT format.
- `collie::BinnedModel` is a ROOT-independent representation for validation and
  likelihood logic.
- `CollieLimit` provides confidence-level calculators and compatibility APIs.
- `collie::ProfileLikelihoodEngine` gives each fit its own Minuit2 state,
  eliminating the former process-global Fortran callback.
- `tools/collie` provides inspection, validation, limits, fits, and combination
  through one interface.

Legacy programs in `examples/` are still built for compatibility; see the
[legacy examples guide](examples/README.md). New code should prefer the unified
CLI and the typed model/run-context APIs.

## Reproducible profile-likelihood example

The self-contained example in `examples/profile_likelihood/` recreates the
four nuisance-profile scenarios and the LLR and CLs scans used in the COLLIE
worked example. It builds and runs in the same pinned ROOT container as CI:

```sh
./examples/profile_likelihood/run.sh
```

The six plots are written as PDF and PNG pairs under
`doc/profile-likelihood-example/`. See the
[example guide](examples/profile_likelihood/README.md) for a quick smoke mode
and configuration options.

The toy model starts with separated background and signal-plus-background NLLR
ensembles, profiles rate and shape nuisance parameters, then scans the signal
strength to obtain the expected and observed CLs limit.

### No nuisance parameters

[![Toy-model background and signal-plus-background NLLR distributions without nuisance parameters][nllr-no-nuisance-png]][nllr-no-nuisance-pdf]

### Fully profiled nuisances

[![Toy-model background and signal-plus-background NLLR distributions with rate and shape nuisances profiled][nllr-profiled-png]][nllr-profiled-pdf]

### LLR scan over signal strength

[![Toy-model LLR scan over the cross-section scale with expected bands][llr-scan-png]][llr-scan-pdf]

### Expected and observed CLs

[![Toy-model expected and observed CLs curves with one- and two-sigma bands][cls-curve-png]][cls-curve-pdf]

[nllr-no-nuisance-png]: ./doc/profile-likelihood-example/profile-likelihood-nllr-no-nuisance-root.png
[nllr-no-nuisance-pdf]: ./doc/profile-likelihood-example/profile-likelihood-nllr-no-nuisance-root.pdf
[nllr-profiled-png]: ./doc/profile-likelihood-example/profile-likelihood-nllr-fully-profiled-root.png
[nllr-profiled-pdf]: ./doc/profile-likelihood-example/profile-likelihood-nllr-fully-profiled-root.pdf
[llr-scan-png]: ./doc/profile-likelihood-example/profile-likelihood-llr-scan-root.png
[llr-scan-pdf]: ./doc/profile-likelihood-example/profile-likelihood-llr-scan-root.pdf
[cls-curve-png]: ./doc/profile-likelihood-example/profile-likelihood-cls-curve-root.png
[cls-curve-pdf]: ./doc/profile-likelihood-example/profile-likelihood-cls-curve-root.pdf
