#!/usr/bin/env bash
# setup.sh -- create .venv, install requirements.txt, and build the
# `dynbem` Rust extension via maturin so `import dynbem` works after
# this script finishes.
#
# Works on Linux / macOS / Windows (git-bash, WSL). Windows callers
# can use the setup.cmd wrapper, which just invokes this script.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"

case "$(uname -s 2>/dev/null || echo unknown)" in
    MINGW*|MSYS*|CYGWIN*) IS_WINDOWS=1 ;;
    *)                     IS_WINDOWS=0 ;;
esac

# ---------------------------------------------------------------------------
# Prerequisite checks. Run them all up front so the user sees every
# missing tool at once instead of failing halfway through.
# ---------------------------------------------------------------------------

PY_REQ_MAJOR=3
PY_REQ_MINOR=10

# Find python (>= 3.10). Prefer `python3` on POSIX, fall back to `python`.
PYTHON_BIN=""
for cand in python3 python; do
    if command -v "$cand" >/dev/null 2>&1; then
        PYTHON_BIN="$cand"
        break
    fi
done

problems=()

if [ -z "$PYTHON_BIN" ]; then
    problems+=("Python interpreter not found on PATH (need python3 or python, >= ${PY_REQ_MAJOR}.${PY_REQ_MINOR}). Install from https://www.python.org/downloads/")
else
    PY_VERSION=$("$PYTHON_BIN" -c 'import sys; print("{}.{}".format(*sys.version_info[:2]))' 2>/dev/null || echo "0.0")
    PY_OK=$("$PYTHON_BIN" -c "import sys; print('1' if sys.version_info >= (${PY_REQ_MAJOR}, ${PY_REQ_MINOR}) else '0')" 2>/dev/null || echo "0")
    if [ "$PY_OK" != "1" ]; then
        problems+=("Python ${PY_VERSION} is too old (need >= ${PY_REQ_MAJOR}.${PY_REQ_MINOR}). Found at: $(command -v "$PYTHON_BIN")")
    fi
fi

if ! command -v cargo >/dev/null 2>&1; then
    problems+=("cargo not on PATH (required to build the dynbem Rust extension via maturin). Install Rust from https://rustup.rs/ and reopen your shell.")
fi

if ! command -v "${CC:-cc}" >/dev/null 2>&1 && ! command -v gcc >/dev/null 2>&1 && ! command -v clang >/dev/null 2>&1 && [ "$IS_WINDOWS" = "0" ]; then
    # POSIX: numpy/numba wheels usually ship binaries but a C compiler
    # is still required to build any source-only deps. On Windows the
    # MSVC toolchain installed alongside Rust handles this.
    problems+=("No C compiler (cc/gcc/clang) on PATH. Install build tools (e.g. 'sudo apt install build-essential' on Debian/Ubuntu, or Xcode CLT on macOS).")
fi

if [ ${#problems[@]} -gt 0 ]; then
    echo "Cannot run setup.sh -- the following prerequisites are missing:" >&2
    echo >&2
    for p in "${problems[@]}"; do
        echo "  * $p" >&2
    done
    echo >&2
    exit 1
fi

echo "Prerequisites:"
echo "  python : $PY_VERSION at $(command -v "$PYTHON_BIN")"
echo "  cargo  : $(cargo --version) at $(command -v cargo)"
echo

# ---------------------------------------------------------------------------
# Step 1: create / find the venv.
# ---------------------------------------------------------------------------

venv_python() {
    if [ -x "$VENV_DIR/bin/python" ]; then
        echo "$VENV_DIR/bin/python"
    elif [ -x "$VENV_DIR/Scripts/python.exe" ]; then
        echo "$VENV_DIR/Scripts/python.exe"
    elif [ -f "$VENV_DIR/Scripts/python.exe" ]; then
        # On some Windows filesystems exec-bit isn't set; -f is enough.
        echo "$VENV_DIR/Scripts/python.exe"
    else
        echo ""
    fi
}

if [ -z "$(venv_python)" ]; then
    echo "Creating virtual environment in $VENV_DIR ..."
    "$PYTHON_BIN" -m venv "$VENV_DIR"
else
    echo "Virtual environment already exists at $VENV_DIR."
fi

VPY="$(venv_python)"
if [ -z "$VPY" ]; then
    echo "error: venv creation failed; no python found under $VENV_DIR" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 2: Python dependencies.
# ---------------------------------------------------------------------------

echo "Upgrading pip ..."
"$VPY" -m pip install --upgrade pip

echo "Installing requirements ..."
"$VPY" -m pip install -r "$SCRIPT_DIR/requirements.txt"

echo "Installing maturin (build tool for the Rust extension) ..."
"$VPY" -m pip install maturin

# ---------------------------------------------------------------------------
# Step 3: build the Rust extension.
# ---------------------------------------------------------------------------

echo "Building the dynbem extension via maturin (release) ..."
# maturin auto-detects the active venv via VIRTUAL_ENV; export it so
# `maturin develop` installs into our .venv even when this script is
# invoked from a parent shell that has its own VIRTUAL_ENV set.
export VIRTUAL_ENV="$VENV_DIR"
(cd "$SCRIPT_DIR/dynbem" && "$VPY" -m maturin develop --release)

echo
echo "Done. Activate the venv with:"
if [ "$IS_WINDOWS" = "1" ]; then
    echo "    .venv\\Scripts\\activate           (cmd / powershell)"
    echo "    source .venv/Scripts/activate    (git-bash)"
else
    echo "    source .venv/bin/activate"
fi
