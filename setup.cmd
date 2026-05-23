@echo off
rem setup.cmd -- thin Windows wrapper around setup.sh.
rem
rem Looks for bash on PATH (e.g. git-bash) and runs setup.sh, which
rem creates the .venv, installs requirements.txt, and builds the
rem dynbem Rust extension via maturin. See setup.sh for details.
setlocal

set "SCRIPT_DIR=%~dp0"

where bash >NUL 2>&1
if errorlevel 1 (
    echo error: bash not found on PATH.
    echo Install Git for Windows ^(https://git-scm.com/download/win^), which
    echo ships git-bash, or use WSL. Then re-run setup.cmd.
    exit /b 1
)

bash "%SCRIPT_DIR%setup.sh" %*
exit /b %ERRORLEVEL%
