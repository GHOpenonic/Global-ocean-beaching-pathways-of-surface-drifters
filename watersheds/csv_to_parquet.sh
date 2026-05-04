#!/bin/bash
set -euo pipefail

INPUT_CSV="${1:-/home/codycruz/drifters_watersheds/undrogued_beach.csv}"
OUTPUT_PARQUET="${2:-/home/codycruz/drifters_watersheds/undrogued_beach.parquet}"
PYTHON_BIN="/home/codycruz/drifters_watersheds/.venv/bin/python"

if [[ ! -f "${INPUT_CSV}" ]]; then
  echo "Input CSV not found: ${INPUT_CSV}" >&2
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "Python environment not found: ${PYTHON_BIN}" >&2
  exit 1
fi

"${PYTHON_BIN}" - <<'PY' "${INPUT_CSV}" "${OUTPUT_PARQUET}"
from pathlib import Path
import sys

input_csv = Path(sys.argv[1])
output_parquet = Path(sys.argv[2])

try:
    import pyarrow.csv as pacsv
    import pyarrow.parquet as papq
except ModuleNotFoundError as exc:
    raise SystemExit(
        "This converter requires pyarrow in /home/codycruz/drifters_watersheds/.venv.\n"
        "Install it with:\n"
        "  /home/codycruz/drifters_watersheds/.venv/bin/pip install pyarrow\n"
        f"Missing module: {exc.name}"
    )

output_parquet.parent.mkdir(parents=True, exist_ok=True)

read_options = pacsv.ReadOptions(block_size=8 * 1024 * 1024)
reader = pacsv.open_csv(input_csv, read_options=read_options)

writer = None
rows_written = 0

try:
    for batch in reader:
        if writer is None:
            writer = papq.ParquetWriter(
                output_parquet,
                batch.schema,
                compression="snappy",
            )
        writer.write_batch(batch)
        rows_written += batch.num_rows
finally:
    if writer is not None:
        writer.close()

if writer is None:
    raise SystemExit(f"No rows were read from {input_csv}")

print(f"Wrote {rows_written} rows to {output_parquet}")
PY
