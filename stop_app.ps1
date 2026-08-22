$ErrorActionPreference = 'Stop'

Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" | Where-Object {
    $_.CommandLine -match 'app\.py'
} | ForEach-Object {
    Stop-Process -Id $_.ProcessId -Force
    Write-Host "Stopped MolSim PID $($_.ProcessId)"
}

if (-not (Get-CimInstance Win32_Process -Filter "Name = 'python.exe'" | Where-Object { $_.CommandLine -match 'app\.py' })) {
    Write-Host "No MolSim process running."
}
