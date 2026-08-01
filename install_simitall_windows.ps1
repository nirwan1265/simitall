param(
    [ValidateSet("minimal", "population", "omics", "sequencing", "full")]
    [string]$Profile = "full",
    [string]$EnvName = "simitall"
)

$ErrorActionPreference = "Stop"
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$Installer = Join-Path $ScriptDir "inst\installers\install_simitall_windows.ps1"
& powershell.exe -NoProfile -ExecutionPolicy Bypass -File $Installer `
    -Profile $Profile -EnvName $EnvName
exit $LASTEXITCODE
