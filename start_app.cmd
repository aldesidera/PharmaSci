@echo off
setlocal
cd /d "%~dp0"
echo [MolSim_ver10] Starting runtime from %CD%
powershell -ExecutionPolicy Bypass -File "%~dp0start_app.ps1"
