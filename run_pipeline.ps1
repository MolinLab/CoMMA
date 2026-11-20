[CmdletBinding()]
param(
    [ValidateRange(1, 1024)] [int]$Workers = 1,
    [ValidateRange(1, 4096)] [double]$MemoryGiB = 2,
    [string]$Config = 'config/breeding_run.yml',
    [string]$RExecutable
)

$ErrorActionPreference = 'Stop'
$root = $PSScriptRoot
if (-not $RExecutable) {
    $install = Get-ItemProperty -LiteralPath 'HKCU:\Software\R-core\R', 'HKLM:\Software\R-core\R' -ErrorAction SilentlyContinue
    $RExecutable = $install | Where-Object { $_.InstallPath } | ForEach-Object { Join-Path $_.InstallPath 'bin\R.exe' } |
        Where-Object { Test-Path -LiteralPath $_ -PathType Leaf } | Select-Object -First 1
}
if (-not $RExecutable) {
    $RExecutable = (Get-Command R.exe -ErrorAction SilentlyContinue | Select-Object -First 1).Source
}
if (-not $RExecutable -or -not (Test-Path -LiteralPath $RExecutable -PathType Leaf)) {
    throw 'Install 64-bit R 4.6.x or pass the path to its R.exe using -RExecutable.'
}
$RExecutable = (Resolve-Path -LiteralPath $RExecutable).Path

Push-Location -LiteralPath $root
try {
    if (-not (Test-Path -LiteralPath $Config -PathType Leaf)) { throw "Configuration file missing: $Config" }
    $env:TEMP = Join-Path $root '.tmp'
    $env:TMP = $env:TEMP
    $env:TMPDIR = $env:TEMP
    $env:R_USER = Join-Path $root '.r_user'
    $env:R_LIBS_USER = Join-Path $root '.breeding-library\R-4.6'
    $env:R_LIBS_SITE = $env:R_LIBS_USER
    $env:R_LIBS = ''
    $env:GWAS_R_EXECUTABLE = $RExecutable
    $env:OMP_NUM_THREADS = '1'
    $env:OPENBLAS_NUM_THREADS = '1'
    $env:MKL_NUM_THREADS = '1'
    $env:BLIS_NUM_THREADS = '1'
    New-Item -ItemType Directory -Force -Path $env:TEMP, $env:R_USER, $env:R_LIBS_USER | Out-Null

    & $RExecutable --vanilla --slave -f 'scripts/setup_breeding_home.R'
    if ($LASTEXITCODE -ne 0) { throw "Dependency setup failed (exit $LASTEXITCODE)." }

    $memory = $MemoryGiB.ToString([System.Globalization.CultureInfo]::InvariantCulture)
    & $RExecutable --vanilla --slave -f 'scripts/run_breeding_home.R' --args $Config $Workers $memory
    if ($LASTEXITCODE -ne 0) { throw "Breeding pipeline failed (exit $LASTEXITCODE)." }
} finally {
    Pop-Location
}
