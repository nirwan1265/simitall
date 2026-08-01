#!/usr/bin/env bash
set -euo pipefail

PROFILE="${1:-full}"
ENV_NAME="${2:-simitall}"
VALID_PROFILES=" minimal population omics sequencing full "
if [[ "$VALID_PROFILES" != *" ${PROFILE} "* ]]; then
  echo "ERROR: profile must be minimal, population, omics, sequencing, or full" >&2
  exit 2
fi

OS="$(uname -s)"
ARCH="$(uname -m)"
if [[ "$OS" != "Darwin" && "$OS" != "Linux" ]]; then
  echo "ERROR: this installer supports macOS, Linux, and Windows through WSL." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/../../DESCRIPTION" ]]; then
  REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
elif [[ -f "${SCRIPT_DIR}/../DESCRIPTION" ]]; then
  REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
else
  echo "ERROR: could not locate the simitall package root." >&2
  exit 2
fi

log() { printf '\n[simitall] %s\n' "$*"; }

find_conda() {
  if command -v conda >/dev/null 2>&1; then
    command -v conda
    return
  fi
  for candidate in \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/mambaforge/bin/conda" \
    "$HOME/miniconda3/bin/conda" \
    "$HOME/anaconda3/bin/conda"; do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return
    fi
  done
  return 1
}

install_miniforge() {
  case "${OS}-${ARCH}" in
    Darwin-arm64) asset="Miniforge3-MacOSX-arm64.sh" ;;
    Darwin-x86_64) asset="Miniforge3-MacOSX-x86_64.sh" ;;
    Linux-x86_64) asset="Miniforge3-Linux-x86_64.sh" ;;
    Linux-aarch64|Linux-arm64) asset="Miniforge3-Linux-aarch64.sh" ;;
    *)
      echo "ERROR: unsupported platform ${OS}-${ARCH}. Install Conda manually." >&2
      exit 2
      ;;
  esac
  local destination="$HOME/miniforge3"
  local installer
  installer="$(mktemp "${TMPDIR:-/tmp}/miniforge.XXXXXX.sh")"
  local url="https://github.com/conda-forge/miniforge/releases/latest/download/${asset}"
  log "Conda was not found; installing Miniforge for ${OS}-${ARCH}."
  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 "$url" -o "$installer"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$installer" "$url"
  else
    echo "ERROR: curl or wget is required to download Miniforge." >&2
    exit 2
  fi
  bash "$installer" -b -p "$destination"
  rm -f "$installer"
  printf '%s\n' "${destination}/bin/conda"
}

CONDA_EXE="$(find_conda || true)"
if [[ -z "$CONDA_EXE" ]]; then
  CONDA_EXE="$(install_miniforge | tail -n 1)"
fi
log "Using Conda: ${CONDA_EXE}"
log "Installation profile: ${PROFILE}"
log "Conda environment: ${ENV_NAME}"

base_packages=(
  "python=3.10" pip "r-base>=4.3" r-remotes r-jsonlite r-matrix
  r-reticulate
)
packages=("${base_packages[@]}")
if [[ "$PROFILE" == "omics" || "$PROFILE" == "full" ]]; then
  packages+=(r-biocmanager r-ggplot2 r-patchwork)
fi
if [[ "$PROFILE" == "sequencing" || "$PROFILE" == "full" ]]; then
  packages+=(art pbsim pbsim3 unicycler samtools seqkit pigz)
fi

if "$CONDA_EXE" env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  log "Updating existing environment ${ENV_NAME}."
  "$CONDA_EXE" install -n "$ENV_NAME" -y --strict-channel-priority \
    -c conda-forge -c bioconda "${packages[@]}"
else
  log "Creating environment ${ENV_NAME}."
  "$CONDA_EXE" create -n "$ENV_NAME" -y --strict-channel-priority \
    -c conda-forge -c bioconda "${packages[@]}"
fi

if [[ "$PROFILE" == "population" || "$PROFILE" == "full" ]]; then
  log "Installing SimuPOP and advanced phenotype dependencies."
  if ! "$CONDA_EXE" install -n "$ENV_NAME" -y -c conda-forge simupop; then
    "$CONDA_EXE" run -n "$ENV_NAME" python -m pip install simuPOP
  fi
  "$CONDA_EXE" run -n "$ENV_NAME" Rscript -e \
    'if (!requireNamespace("simplePHENOTYPES", quietly=TRUE)) remotes::install_github("samuelbfernandes/simplePHENOTYPES", dependencies=TRUE, upgrade="never")'
fi

if [[ "$PROFILE" == "omics" || "$PROFILE" == "full" ]]; then
  log "Installing Bioconductor omics packages."
  "$CONDA_EXE" run -n "$ENV_NAME" Rscript -e \
    'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("Rsubread", "ChIPsim", "SingleCellExperiment", "SummarizedExperiment", "splatter"), ask=FALSE, update=FALSE)'
fi

if [[ "$PROFILE" == "sequencing" || "$PROFILE" == "full" ]]; then
  log "Installing Python sequencing tools."
  "$CONDA_EXE" run -n "$ENV_NAME" python -m pip install badread quast
fi

log "Installing the simitall R package."
export SIMITALL_REPO_ROOT="$REPO_ROOT"
"$CONDA_EXE" run -n "$ENV_NAME" Rscript -e \
  'remotes::install_local(Sys.getenv("SIMITALL_REPO_ROOT"), dependencies=FALSE, upgrade="never", force=TRUE)'

log "Checking installed dependencies."
"$CONDA_EXE" run -n "$ENV_NAME" Rscript -e \
  "simitall::check_simitall_dependencies('${PROFILE}')"

cat <<EOF

[simitall] Installation complete.

Activate the environment with:
  conda activate ${ENV_NAME}

Then start R and run:
  library(simitall)
  check_simitall_dependencies("${PROFILE}")
EOF
