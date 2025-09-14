#!/usr/bin/env bash
# Create a clean submission zip with code + build scripts, no builds/executables.
# Usage: scripts/make_submission.sh [output.zip]
set -euo pipefail

# Absolute path to this project root
ROOT="/home/sghys/projects/CPSC524-modeling"
cd "$ROOT"

OUT="${1:-submission_$(date +%Y%m%d_%H%M%S).zip}"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

# ---- WHAT TO INCLUDE (edit if you need more) ----
INCLUDES=(
  CMakeLists.txt
  cmake/                   # your cmake helper macros
  src/                     # all source (library, apps, tests if you keep them)
  third-party/doctest/     # the single-header doctest (if you vendored it)
  README*                  # README file(s) with your name & student number
)

# ---- WHAT TO EXCLUDE (assignment rules + common junk) ----
# You can tweak this list, but it already covers the assignment’s requirements.
read -r -d '' RSYNC_EXCLUDES <<'EOF' || true
--exclude .git/
--exclude .gitmodules
--exclude .github/
--exclude .idea/
--exclude .vs/
--exclude build*/
--exclude CMakeFiles/
--exclude CMakeCache.txt
--exclude Testing/
--exclude _deps/
--exclude compile_commands.json
--exclude **/*.o
--exclude **/*.obj
--exclude **/*.a
--exclude **/*.lib
--exclude **/*.so
--exclude **/*.dylib
--exclude **/*.dll
--exclude **/*.exe
--exclude **/*.pdb
--exclude **/*.pyc
--exclude __pycache__/
--exclude .DS_Store
EOF

# copy includes into staging with excludes applied
for path in "${INCLUDES[@]}"; do
  if [ -e "$path" ]; then
    rsync -a ${RSYNC_EXCLUDES} "$path" "$STAGE"/
  fi
done

echo "[i] PWD:  $(pwd)"
echo "[i] ROOT: ${ROOT:-<unset>}"
echo "[i] STAGE:${STAGE:-<unset>}"

echo "[i] INCLUDES (verbatim):"

echo "[i] Files staged:"
find "$STAGE" -type f | sed "s#${STAGE}#.#"

# build the zip (needs 'zip' installed in WSL: sudo apt-get install zip)
cd "$STAGE"
zip -r -9 "$OUT" . >/dev/null
mv "$OUT" "$ROOT"/

echo "Created: $ROOT/$OUT"
