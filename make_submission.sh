#!/usr/bin/env bash
# Create a clean submission zip with code + build scripts, no builds/executables.
# Always writes ZIP to ~/projects/CPSC524-modeling
set -euo pipefail

ROOT="${HOME}/projects/CPSC524-modeling"
if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: Root not found: $ROOT" >&2
  exit 1
fi
cd "$ROOT"

# Output goes to ROOT, name can be overridden by arg (basename only)
BASENAME="${1:-submission_$(date +%Y%m%d_%H%M%S).zip}"
OUT_PATH="$ROOT/$(basename "$BASENAME")"

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# ----- INCLUDE EXACTLY WHAT WE WANT -----
INCLUDES=(
  CMakeLists.txt
  cmake               # CMake modules/macros
  include             # (if present)
  src                 # ALL source
  mesh                # ALL meshes
  scripts             # helper scripts (including this one)
  third-party         # (if allowed; remove if not)
  README*             # README files
)

# ----- EXCLUDES (no builds/binaries/IDE junk) -----
RSYNC_EXCLUDES=(
  --exclude .git/
  --exclude .gitmodules
  --exclude .github/
  --exclude .idea/
  --exclude .vscode/
  --exclude .vs/
  --exclude 'build*/'
  --exclude CMakeFiles/
  --exclude CMakeCache.txt
  --exclude Testing/
  --exclude _deps/
  --exclude compile_commands.json
  # object files & libs
  --exclude '**/*.o'
  --exclude '**/*.a'
  --exclude '**/*.so'
  --exclude '**/*.dylib'
  --exclude '**/*.dll'
  --exclude '**/*.lib'
  # Windows/MSVC .obj ONLY in build artifacts (keep mesh/*.obj)
  --exclude 'build*/**/*.obj'
  --exclude 'CMakeFiles/**/*.obj'
  # executables & debug
  --exclude '**/*.exe'
  --exclude '**/*.pdb'
  # py/OS cruft
  --exclude '**/*.pyc'
  --exclude __pycache__/
  --exclude .DS_Store
)

# Stage files
for path in "${INCLUDES[@]}"; do
  if [[ -e "$path" ]]; then
    rsync -a --prune-empty-dirs "${RSYNC_EXCLUDES[@]}" "$path" "$STAGE"/
  fi
done

# Must have full src/ and mesh/
if [[ ! -d "$STAGE/src" ]]; then
  echo "ERROR: src/ not staged. Check path: $ROOT/src" >&2
  exit 2
fi
if [[ ! -d "$STAGE/mesh" ]]; then
  echo "ERROR: mesh/ not staged. Check path: $ROOT/mesh" >&2
  exit 3
fi

# Ensure a README exists (create a template if missing)
if ! find "$STAGE" -maxdepth 1 -iname 'README*' | grep -q .; then
  cat > "$STAGE/README_Submission.txt" <<'TXT'
# Assignment Submission

Name: <YOUR NAME HERE>
Student Number: <YOUR STUDENT NUMBER HERE>

Build Instructions:
- mkdir -p build && cd build
- cmake ..
- cmake --build .
TXT
fi

echo "[i] Staging complete. File count:"
( cd "$STAGE" && find . -type f | wc -l )

# Build zip into ROOT
( cd "$STAGE" && zip -r -9 "$(basename "$OUT_PATH")" . >/dev/null )
mv "$STAGE/$(basename "$OUT_PATH")" "$OUT_PATH"

echo "Created: $OUT_PATH"
