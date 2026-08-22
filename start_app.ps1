$ErrorActionPreference = 'Stop'

$projectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$logsDir = Join-Path $projectDir 'logs'
If (-not (Test-Path $logsDir)) {
    New-Item -ItemType Directory -Path $logsDir -Force | Out-Null
}

$existing = Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" | Where-Object {
    $_.CommandLine -match 'app\.py'
}

If ($existing) {
    Write-Host "MolSim already running on PID $($existing[0].ProcessId)"
    exit 0
}

$logOut = Join-Path $logsDir 'app.out.log'
$logErr = Join-Path $logsDir 'app.err.log'

Start-Process -FilePath 'python' -ArgumentList 'app.py' -WorkingDirectory $projectDir -WindowStyle Hidden -RedirectStandardOutput $logOut -RedirectStandardError $logErr | Out-Null

Write-Host "MolSim started in $projectDir"
Write-Host "Logs: $logOut / $logErr"
