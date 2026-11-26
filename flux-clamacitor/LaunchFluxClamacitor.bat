@echo off
title Flux Clamacitor Launcher
echo.
echo  ======================================================
echo   FLUX CLAMACITOR
echo   CaCO3 Removal Flux to Population Estimator
echo   Corbicula fluminea Mass Balance Model
echo  ======================================================
echo.
echo  Launching application...
echo.

:: Get the directory this .bat file lives in
set "SCRIPT_DIR=%~dp0"
set "HTML_FILE=%SCRIPT_DIR%FluxClamacitor.html"

:: Check the HTML file exists
if not exist "%HTML_FILE%" (
    echo  ERROR: FluxClamacitor.html not found in the same folder as this .bat file.
    echo  Expected: %HTML_FILE%
    pause
    exit /b 1
)

:: Try to open with the default browser
start "" "%HTML_FILE%"

echo  Browser launched. If nothing opened, open FluxClamacitor.html manually.
echo.
timeout /t 3 /nobreak >nul
exit /b 0
