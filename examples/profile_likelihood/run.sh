#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
image="${ROOT_IMAGE:-rootproject/root:6.30.06-ubuntu22.04}"
platform="${ROOT_PLATFORM:-linux/amd64}"
effort="${COLLIE_EFFORT:-veryfine}"
build_jobs="${COLLIE_BUILD_JOBS:-4}"
output_dir="${COLLIE_OUTPUT_DIR:-$repo_root/doc/profile-likelihood-example}"

case "$effort" in
  test|veryveryfast|veryfast|fast|standard|fine|veryfine) ;;
  *)
    echo "Unsupported COLLIE_EFFORT: $effort" >&2
    exit 1
    ;;
esac

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is required to run the pinned ROOT environment" >&2
  exit 1
fi

mkdir -p "$output_dir"
output_dir="$(cd -- "$output_dir" && pwd)"

remove_work_dir=false
if [[ -n "${COLLIE_WORK_DIR:-}" ]]; then
  mkdir -p "$COLLIE_WORK_DIR"
  work_dir="$(cd -- "$COLLIE_WORK_DIR" && pwd)"
else
  work_dir="$(mktemp -d /tmp/collie-profile-likelihood.XXXXXX)"
  remove_work_dir=true
fi

echo "Building the profile-likelihood example with $image"
echo "Scenario effort: $effort"
echo "Plot destination: $output_dir"
echo "Restart state: $work_dir"

if ! docker run --rm \
  --platform "$platform" \
  -e "BUILD_JOBS=$build_jobs" \
  -e "COLLIE_EFFORT=$effort" \
  -e "PROFILE_LIKELIHOOD_COLLIE_RESULTS=/state/results-$effort" \
  -e PROFILE_LIKELIHOOD_COLLIE_OUTPUT=/output \
  -v "$repo_root:/source:ro" \
  -v "$output_dir:/output" \
  -v "$work_dir:/state" \
  -w /source \
  "$image" \
  bash -lc '
    set -euo pipefail

    build_dir=/state/build
    build_log=/state/build.log
    if ! cmake -S /source -B "$build_dir" \
      -DCMAKE_BUILD_TYPE=Release \
      -DCOLLIE_BUILD_EXAMPLES=ON \
      -DBUILD_TESTING=OFF >"$build_log" 2>&1; then
      tail -n 120 "$build_log" >&2
      exit 1
    fi
    if ! cmake --build "$build_dir" \
      --target collieProfileLikelihood \
      --parallel "$BUILD_JOBS" >>"$build_log" 2>&1; then
      tail -n 120 "$build_log" >&2
      exit 1
    fi
    echo "Built collieProfileLikelihood"

    driver="$build_dir/bin/collieProfileLikelihood.exe"
    result_dir="$PROFILE_LIKELIHOOD_COLLIE_RESULTS"
    mkdir -p "$result_dir"
    "$driver" create "$result_dir/profile-likelihood-toy-input.root"

    seeds=(1747 1781 1811 1847)
    pids=()
    scenario_start=$SECONDS
    for scenario in 0 1 2 3; do
      scenario_result="$result_dir/scenario-$scenario.root"
      scenario_done="$result_dir/scenario-$scenario.done"
      if [[ -s "$scenario_result" && -f "$scenario_done" ]]; then
        echo "Reusing completed scenario $scenario"
        continue
      fi
      (
        scenario_work="/tmp/collie-scenario-$scenario"
        mkdir -p "$scenario_work"
        cd "$scenario_work"
        "$driver" \
          scenario \
          "$result_dir/profile-likelihood-toy-input.root" \
          "$scenario_result" \
          "$scenario" \
          "$COLLIE_EFFORT" \
          "${seeds[$scenario]}"
        touch "$scenario_done"
      ) >"$result_dir/scenario-$scenario.log" 2>&1 &
      pids+=("$!")
    done

    next_report=30
    while :; do
      running=0
      for pid in "${pids[@]}"; do
        if kill -0 "$pid" 2>/dev/null; then
          running=1
        fi
      done
      if [[ "$running" -eq 0 ]]; then
        break
      fi
      sleep 5
      elapsed=$((SECONDS - scenario_start))
      if [[ "$elapsed" -ge "$next_report" ]]; then
        echo "COLLIE profile ensembles are still running (${elapsed}s)..."
        next_report=$((next_report + 30))
      fi
    done

    status=0
    for pid in "${pids[@]}"; do
      if ! wait "$pid"; then
        status=1
      fi
    done
    if [[ "$status" -ne 0 ]]; then
      for log in "$result_dir"/scenario-*.log; do
        echo "===== $log =====" >&2
        tail -n 80 "$log" >&2
      done
      exit 1
    fi

    grep -h "Finished" "$result_dir"/scenario-*.log || true
    fit_warnings=$({
      grep -h "Fit did not converge\|Fit didn.t converge" \
        "$result_dir"/scenario-*.log || true
    } | wc -l | tr -d " ")
    echo "Rejected COLLIE fit attempts: $fit_warnings"
    echo "COLLIE profile ensembles completed in $((SECONDS - scenario_start))s"

    echo "Running the COLLIE AWW scan..."
    if [[ -s "$result_dir/scan.root" && -f "$result_dir/scan.done" ]]; then
      echo "Reusing completed COLLIE scan"
    else
      if ! "$driver" \
        scan \
        "$result_dir/profile-likelihood-toy-input.root" \
        "$result_dir/scan.root" >"$result_dir/scan.log" 2>&1; then
        tail -n 120 "$result_dir/scan.log" >&2
        exit 1
      fi
      touch "$result_dir/scan.done"
    fi
    if grep -Eiq \
      "Fit did.?n.t converge|(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|barf" \
      "$result_dir/scan.log"; then
      echo "COLLIE AWW scan reported an invalid fit" >&2
      grep -Ein \
        "Fit did.?n.t converge|(^|[^[:alpha:]])(nan|inf)([^[:alpha:]]|$)|barf" \
        "$result_dir/scan.log" >&2
      exit 1
    fi
    grep "Finished" "$result_dir/scan.log" || true

    root -l -b -q \
      /source/examples/profile_likelihood/profile_likelihood_collie_plots.C

    stems=(
      profile-likelihood-nllr-no-nuisance-root
      profile-likelihood-nllr-rate-constrained-root
      profile-likelihood-nllr-shape-constrained-root
      profile-likelihood-nllr-fully-profiled-root
      profile-likelihood-llr-scan-root
      profile-likelihood-cls-curve-root
    )
    for stem in "${stems[@]}"; do
      test -s "/output/$stem.pdf"
      test -s "/output/$stem.png"
    done
  '; then
  echo "Run failed; restart state retained at $work_dir" >&2
  exit 1
fi

find "$output_dir" -maxdepth 1 -type f \
  \( -name 'profile-likelihood-*.pdf' -o -name 'profile-likelihood-*.png' \) \
  -print | sort

if [[ "$remove_work_dir" == true \
  && "$work_dir" == /tmp/collie-profile-likelihood.* ]]; then
  find "$work_dir" -depth -delete
fi
