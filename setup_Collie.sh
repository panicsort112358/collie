SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  # Resolve symlinks until SOURCE names the underlying file.
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  if [[ $SOURCE != /* ]]; then
    SOURCE="$DIR/$SOURCE"
  fi
done
COLLIE="$(cd -P "$(dirname "$SOURCE")" && pwd)"

if ! command -v root-config >/dev/null 2>&1; then
  echo "ERROR: root-config is not on PATH; source your ROOT environment first" >&2
  return 1 2>/dev/null || exit 1
fi

export LD_LIBRARY_PATH="$COLLIE/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$COLLIE/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
export PATH="$COLLIE/tools:$COLLIE/examples:$PATH"
