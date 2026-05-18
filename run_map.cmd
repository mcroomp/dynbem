@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "VENV_PY=%SCRIPT_DIR%.venv\Scripts\python.exe"

if not exist "%VENV_PY%" (
    echo Virtual environment not found at %VENV_PY%
    echo Run setup.cmd first.
    exit /b 1
)

REM Default: quick grid, save to out\map.npz, plot to out\
REM Pass any args (e.g. --full, --quantity rpm, --tmax 800) to override.

if "%~1"=="" (
    if not exist "%SCRIPT_DIR%out" mkdir "%SCRIPT_DIR%out"
    "%VENV_PY%" -m envelope.compute_map --quick --save "%SCRIPT_DIR%out\map.npz" --plot "%SCRIPT_DIR%out"
) else (
    "%VENV_PY%" -m envelope.compute_map %*
)

endlocal
