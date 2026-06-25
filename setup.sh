#!/usr/bin/env bash
# setup.sh -- bootstrap the dev environment with uv.
#
# uv owns the environment for this repo: it creates the .venv, resolves
# the uv + Cargo workspaces, and builds the `dynbem` Rust extension via
# maturin so `import dynbem` works after this script finishes.
#
# Works on Linux / macOS / Windows (git-bash, WSL). Windows callers can
# use the setup.cmd wrapper, which runs `uv sync` directly (no bash
# required).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Prerequisite checks. Run them all up front so the user sees every
# missing tool at once instead of failing halfway through.
# ---------------------------------------------------------------------------

problems=()

if ! command -v uv >/dev/null 2>&1; then
    problems+=("uv not on PATH. Install it with 'curl -LsSf https://astral.sh/uv/install.sh | sh' (POSIX) or see https://docs.astral.sh/uv/getting-started/installation/ , then reopen your shell.")
fi

if ! command -v cargo >/dev/null 2>&1; then
    problems+=("cargo not on PATH (required to build the dynbem Rust extension via maturin). Install Rust from https://rustup.rs/ and reopen your shell.")
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
echo "  uv    : $(uv --version) at $(command -v uv)"
echo "  cargo : $(cargo --version) at $(command -v cargo)"
echo

# ---------------------------------------------------------------------------
# Sync the workspace. `uv sync --group dev` creates .venv (if needed),
# installs runtime + dev dependencies, and builds the dynbem extension
# editable via maturin against the sibling dynbem_rs crate.
# ---------------------------------------------------------------------------

echo "Syncing environment (uv sync --group dev) ..."
(cd "$SCRIPT_DIR" && uv sync --group dev)

echo
echo "Done. Run commands through uv, e.g.:"
echo "    uv run pytest tests/ -q"
echo "    uv run python -m envelope.compute_map --help"
echo
echo "Or activate the uv-managed venv directly:"
case "$(uname -s 2>/dev/null || echo unknown)" in
    MINGW*|MSYS*|CYGWIN*)
        echo "    .venv\\Scripts\\activate           (cmd / powershell)"
        echo "    source .venv/Scripts/activate    (git-bash)"
        ;;
    *)
        echo "    source .venv/bin/activate"
        ;;
esac
