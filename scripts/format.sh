#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
mode="${1:-format}"

case "$mode" in
  format)
    clang_format_args=(-i)
    ;;
  --check)
    clang_format_args=(--dry-run --Werror)
    ;;
  *)
    echo "usage: $0 [--check]" >&2
    exit 2
    ;;
esac

if ! command -v clang-format >/dev/null 2>&1; then
  echo "clang-format is required" >&2
  exit 1
fi

while IFS= read -r -d '' source; do
  clang-format \
    --style=file \
    --fail-on-incomplete-format \
    "${clang_format_args[@]}" \
    "$repo_root/$source"
done < <(
  git -C "$repo_root" ls-files -z --cached --others --exclude-standard -- \
    '*.C' '*.cc' '*.cpp' '*.cxx' '*.h' '*.hh' '*.hpp'
)
