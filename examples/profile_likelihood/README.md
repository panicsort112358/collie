# Profile-likelihood worked example

This example reproduces the COLLIE plots used by the profile-likelihood toy in
the doctoral-thesis application. It is fully local to this repository: the
25-bin input model is generated analytically, built against the current COLLIE
libraries, and rendered with ROOT without reading the thesis checkout.

Run the full reproduction from the repository root:

```sh
./examples/profile_likelihood/run.sh
```

The default `veryfine` effort generates 200,000 pseudo-experiments per
hypothesis for each of four scenarios. For a quick end-to-end smoke run, use:

```sh
COLLIE_EFFORT=test ./examples/profile_likelihood/run.sh
```

Runtime is dominated by the serial `CLfit2` ensembles. On the tested
Apple-silicon Docker setup, the full ensemble phase took 3,403 seconds
(56 minutes 43 seconds); the smoke mode completed the same phase in 5 seconds.
Cold container builds add roughly one minute, and timings will vary by host.

The runner uses `rootproject/root:6.30.06-ubuntu22.04` on `linux/amd64`, builds
in a temporary container directory, and writes only the finished figures to
`doc/profile-likelihood-example/`. Override its settings with:

- `ROOT_IMAGE` and `ROOT_PLATFORM` for the container runtime.
- `COLLIE_BUILD_JOBS` for build parallelism.
- `COLLIE_EFFORT` for `test`, `veryveryfast`, `veryfast`, `fast`, `standard`,
  `fine`, or `veryfine`.
- `COLLIE_OUTPUT_DIR` for a different host output directory.
- `COLLIE_WORK_DIR` for persistent build/results state that can resume an
  interrupted run.

Each run creates these PDF/PNG pairs:

- `profile-likelihood-nllr-no-nuisance-root`
- `profile-likelihood-nllr-rate-constrained-root`
- `profile-likelihood-nllr-shape-constrained-root`
- `profile-likelihood-nllr-fully-profiled-root`
- `profile-likelihood-llr-scan-root`
- `profile-likelihood-cls-curve-root`

The four ensembles use fixed seeds and the scan uses a fixed seed, so repeated
runs in the same standard-library environment are reproducible. By default,
intermediate ROOT files and logs are removed after a successful run; setting
`COLLIE_WORK_DIR` preserves them for restart and inspection.
