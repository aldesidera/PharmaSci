$ErrorActionPreference = 'Stop'

$projectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectName = Split-Path $projectDir -Leaf
$mainScript = Join-Path $projectDir 'main.py'

$running = Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" | Where-Object {
    $_.CommandLine -match [Regex]::Escape($mainScript)
}

if ($running) {
    $running | ForEach-Object {
        Stop-Process -Id $_.ProcessId -Force
        Write-Host "[$projectName] Stopped runtime PID $($_.ProcessId)"
    }
    exit 0
}

Write-Host "[$projectName] No runtime process running in $projectDir."
