@echo off
rem setup.cmd -- Windows bootstrap for the dev environment using uv.
rem
rem uv owns the environment: it creates the .venv, resolves the uv +
rem Cargo workspaces, and builds the `dynbem` Rust extension via maturin.
rem This wrapper runs `uv sync` directly -- no bash/git-bash required.
setlocal

set "SCRIPT_DIR=%~dp0"

where uv >NUL 2>&1
if errorlevel 1 (
    echo error: uv not found on PATH.
    echo Install it with:  powershell -c "irm https://astral.sh/uv/install.ps1 ^| iex"
    echo or see https://docs.astral.sh/uv/getting-started/installation/
    echo Then reopen your shell and re-run setup.cmd.
    exit /b 1
)

where cargo >NUL 2>&1
if errorlevel 1 (
    echo error: cargo not found on PATH.
    echo The dynbem Rust extension is built via maturin, which needs a Rust
    echo toolchain. Install Rust from https://rustup.rs/ and reopen your shell.
    exit /b 1
)

echo Syncing environment (uv sync --group dev) ...
pushd "%SCRIPT_DIR%"
uv sync --group dev
set "RC=%ERRORLEVEL%"
popd
if not "%RC%"=="0" exit /b %RC%

echo.
echo Done. Run commands through uv, e.g.:
echo     uv run pytest tests/ -q
echo     uv run python -m envelope.compute_map --help
exit /b 0
