@echo off
setlocal EnableExtensions

rem Build a Sigma17{223}[110] FCC Cu symmetric tilt grain boundary with Atomsk.
rem Grain 2 uses the Al label only to retain LAMMPS atom type 2. The LAMMPS
rem input maps both types to a Cu EAM potential during energy minimization.

set "ATOMSK_BIN=D:\2025.916aftersoftware\atmosk\Atomsk\atomsk.exe"
set "LATTICE_CONSTANT=3.597"
set "DUPLICATE_X=20"
set "DUPLICATE_Y=20"
set "DUPLICATE_Z=1"
set "OUTPUT=sigma17_stgb_initial.lmp"

if not exist "%ATOMSK_BIN%" (
  echo ERROR: Atomsk executable was not found: %ATOMSK_BIN%
  echo Edit ATOMSK_BIN in this file to match your local Atomsk installation.
  exit /b 1
)

"%ATOMSK_BIN%" --create fcc %LATTICE_CONSTANT% Cu orient [-33-4] [2-2-3] [110] -duplicate %DUPLICATE_X% %DUPLICATE_Y% %DUPLICATE_Z% grain_1.xsf
if errorlevel 1 exit /b 1

"%ATOMSK_BIN%" --create fcc %LATTICE_CONSTANT% Al orient [-334] [-22-3] [110] -duplicate %DUPLICATE_X% %DUPLICATE_Y% %DUPLICATE_Z% grain_2.xsf
if errorlevel 1 exit /b 1

"%ATOMSK_BIN%" --merge Y 2 grain_1.xsf grain_2.xsf %OUTPUT%
if errorlevel 1 exit /b 1

echo Created %OUTPUT%
echo Copy %OUTPUT% into the lammps folder before running the LAMMPS inputs.
