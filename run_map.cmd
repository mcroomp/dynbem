@echo off
setlocal

set "SCRIPT_DIR=%~dp0"

where uv >NUL 2>&1
if errorlevel 1 (
    echo error: uv not found on PATH. Run setup.cmd first.
    exit /b 1
)

REM Default: quick grid, save to out\map.npz, plot to out\
REM Pass any args (e.g. --full, --quantity rpm, --tmax 800) to override.

pushd "%SCRIPT_DIR%"
if "%~1"=="" (
    if not exist "%SCRIPT_DIR%out" mkdir "%SCRIPT_DIR%out"
    uv run python -m envelope.compute_map --quick --save "%SCRIPT_DIR%out\map.npz" --plot "%SCRIPT_DIR%out"
) else (
    uv run python -m envelope.compute_map %*
)
set "RC=%ERRORLEVEL%"
popd

endlocal & exit /b %RC%
