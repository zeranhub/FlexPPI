@echo off
setlocal
set "CASE_PYMOL_RUNTIME=D:\Anaconda\envs\pvb_cpu"
if not exist "%CASE_PYMOL_RUNTIME%\pythonw.exe" (
    echo Please open ubiquitin/ubiquitin_2QHO.pse with your installed PyMOL.
    pause
    exit /b 1
)
set "PATH=%CASE_PYMOL_RUNTIME%;%CASE_PYMOL_RUNTIME%\Library\bin;%CASE_PYMOL_RUNTIME%\Scripts;%PATH%"
start "" "%CASE_PYMOL_RUNTIME%\pythonw.exe" -m pymol "%~dp0ubiquitin\ubiquitin_2QHO.pse"
endlocal
