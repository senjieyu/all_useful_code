@echo off
setlocal
cd /d "%~dp0"

REM Package the double-clickable launcher. Keep configs external so parameters
REM can be edited without rebuilding the executable.
pyinstaller --noconfirm --clean --onefile --windowed --name CSLGrainBoundaryBuilder --paths src run_builder.py
if errorlevel 1 exit /b %errorlevel%

if not exist "dist\configs" mkdir "dist\configs"
copy /Y "configs\sigma17_110_223_cu.json" "dist\configs\sigma17_110_223_cu.json" >nul
echo.
echo Done: dist\CSLGrainBoundaryBuilder.exe
echo Edit: dist\configs\sigma17_110_223_cu.json
