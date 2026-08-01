param(
    [ValidateSet("minimal", "population", "omics", "sequencing", "full")]
    [string]$Profile = "full",
    [string]$EnvName = "simitall"
)

$ErrorActionPreference = "Stop"

function Write-Step([string]$Message) {
    Write-Host "`n[simitall] $Message" -ForegroundColor Cyan
}

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$SourceCandidate = Resolve-Path (Join-Path $ScriptDir "..\..")
if (Test-Path (Join-Path $SourceCandidate "DESCRIPTION")) {
    $RepoRoot = $SourceCandidate.Path
} else {
    $InstalledCandidate = Resolve-Path (Join-Path $ScriptDir "..")
    if (-not (Test-Path (Join-Path $InstalledCandidate "DESCRIPTION"))) {
        throw "Could not locate the simitall package root."
    }
    $RepoRoot = $InstalledCandidate.Path
}

# ART, PBSIM, and Unicycler are distributed through Bioconda, which does not
# provide native Windows packages. Delegate these profiles to Linux in WSL2.
if ($Profile -in @("sequencing", "full")) {
    Write-Step "The $Profile profile uses WSL2 for Linux/Bioconda sequencing tools."
    $Wsl = Get-Command wsl.exe -ErrorAction SilentlyContinue
    if (-not $Wsl) {
        throw "WSL2 is required for '$Profile'. Run 'wsl --install -d Ubuntu', restart Windows, and run this installer again."
    }
    $WslPath = (& wsl.exe wslpath -a $RepoRoot).Trim()
    if ($LASTEXITCODE -ne 0 -or -not $WslPath) {
        throw "WSL is installed but no Linux distribution is ready. Install Ubuntu with 'wsl --install -d Ubuntu'."
    }
    & wsl.exe bash -lc "cd '$WslPath' && bash ./install_simitall_linux.sh '$Profile' '$EnvName'"
    exit $LASTEXITCODE
}

function Find-Conda {
    $Command = Get-Command conda.exe -ErrorAction SilentlyContinue
    if ($Command) { return $Command.Source }
    $Candidates = @(
        "$env:USERPROFILE\miniforge3\Scripts\conda.exe",
        "$env:USERPROFILE\miniconda3\Scripts\conda.exe",
        "$env:USERPROFILE\anaconda3\Scripts\conda.exe",
        "$env:LOCALAPPDATA\miniforge3\Scripts\conda.exe"
    )
    foreach ($Candidate in $Candidates) {
        if (Test-Path $Candidate) { return $Candidate }
    }
    return $null
}

$Conda = Find-Conda
if (-not $Conda) {
    Write-Step "Conda was not found; installing Miniforge for Windows x86-64."
    $Installer = Join-Path $env:TEMP "Miniforge3-Windows-x86_64.exe"
    $Url = "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe"
    Invoke-WebRequest -Uri $Url -OutFile $Installer
    $Destination = Join-Path $env:USERPROFILE "miniforge3"
    Start-Process -FilePath $Installer -ArgumentList @("/S", "/InstallationType=JustMe", "/AddToPath=0", "/D=$Destination") -Wait
    Remove-Item $Installer -Force
    $Conda = Join-Path $Destination "Scripts\conda.exe"
}
if (-not (Test-Path $Conda)) { throw "Conda installation was not found." }

Write-Step "Using Conda: $Conda"
Write-Step "Installation profile: $Profile"

$Packages = @(
    "python=3.10", "pip", "r-base>=4.3", "r-remotes", "r-jsonlite",
    "r-matrix", "r-reticulate"
)
if ($Profile -eq "omics") {
    $Packages += @("r-biocmanager", "r-ggplot2", "r-patchwork")
}

$Environments = & $Conda env list
if ($Environments -match "(?m)^$([regex]::Escape($EnvName))\s") {
    Write-Step "Updating existing environment $EnvName."
    & $Conda install -n $EnvName -y --strict-channel-priority -c conda-forge $Packages
} else {
    Write-Step "Creating environment $EnvName."
    & $Conda create -n $EnvName -y --strict-channel-priority -c conda-forge $Packages
}
if ($LASTEXITCODE -ne 0) { throw "Conda environment creation failed." }

if ($Profile -eq "population") {
    Write-Step "Installing SimuPOP and advanced phenotype dependencies."
    & $Conda install -n $EnvName -y -c conda-forge simupop
    if ($LASTEXITCODE -ne 0) {
        & $Conda run -n $EnvName python -m pip install simuPOP
    }
    & $Conda run -n $EnvName Rscript -e 'if (!requireNamespace("simplePHENOTYPES", quietly=TRUE)) remotes::install_github("samuelbfernandes/simplePHENOTYPES", dependencies=TRUE, upgrade="never")'
}

if ($Profile -eq "omics") {
    Write-Step "Installing Bioconductor omics packages."
    & $Conda run -n $EnvName Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("Rsubread", "ChIPsim", "SingleCellExperiment", "SummarizedExperiment", "splatter"), ask=FALSE, update=FALSE)'
}

Write-Step "Installing the simitall R package."
$env:SIMITALL_REPO_ROOT = $RepoRoot
& $Conda run -n $EnvName Rscript -e 'remotes::install_local(Sys.getenv("SIMITALL_REPO_ROOT"), dependencies=FALSE, upgrade="never", force=TRUE)'
if ($LASTEXITCODE -ne 0) { throw "The simitall R package installation failed." }

Write-Step "Checking installed dependencies."
& $Conda run -n $EnvName Rscript -e "simitall::check_simitall_dependencies('$Profile')"

Write-Host "`n[simitall] Installation complete." -ForegroundColor Green
Write-Host "Open a Miniforge Prompt and run: conda activate $EnvName"
