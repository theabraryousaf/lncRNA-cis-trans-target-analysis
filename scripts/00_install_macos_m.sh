#!/usr/bin/env bash
set -euo pipefail

# Install command-line tools on macOS Apple Silicon.
brew install micromamba git wget curl bedtools

# Initialize micromamba for zsh if needed.
# Run only once if micromamba activation does not work yet.
micromamba shell init -s zsh -r ~/micromamba || true

echo "Now run: source ~/.zshrc"
echo "Then create the environment:"
echo "micromamba create -y -n lncrna_cistrans -f environment-macos-m.yml"
echo "micromamba activate lncrna_cistrans"
