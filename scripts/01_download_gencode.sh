#!/usr/bin/env bash
set -euo pipefail

source config/config.env
mkdir -p ref

if [[ -s "$GTF_FILE" ]]; then
  echo "GTF already exists: $GTF_FILE"
  exit 0
fi

GZ_FILE="${GTF_FILE}.gz"

echo "Downloading GENCODE GTF..."
curl -L -o "$GZ_FILE" "$GENCODE_GTF_URL"

echo "Uncompressing GTF..."
gunzip -c "$GZ_FILE" > "$GTF_FILE"

ls -lh "$GTF_FILE"
