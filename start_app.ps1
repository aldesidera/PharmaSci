$ErrorActionPreference = 'Stop'

$projectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectName = Split-Path $projectDir -Leaf
$mainScript = Join-Path $projectDir 'main.py'
$logsDir = Join-Path $projectDir 'logs'
If (-not (Test-Path $logsDir)) {
    New-Item -ItemType Directory -Path $logsDir -Force | Out-Null
}

$existing = Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" | Where-Object {
    $_.CommandLine -match [Regex]::Escape($mainScript)
}

If ($existing) {
    Write-Host "[$projectName] Runtime already running on PID $($existing[0].ProcessId)"
    exit 0
}

$logOut = Join-Path $logsDir 'app.out.log'
$logErr = Join-Path $logsDir 'app.err.log'

Start-Process -FilePath 'python' -ArgumentList $mainScript -WorkingDirectory $projectDir -WindowStyle Hidden -RedirectStandardOutput $logOut -RedirectStandardError $logErr | Out-Null

Write-Host "[$projectName] Runtime ativo: $projectDir"
Write-Host "[$projectName] Start URL: http://127.0.0.1:5000"
Write-Host "[$projectName] Logs: $logOut / $logErr"
