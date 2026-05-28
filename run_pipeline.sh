#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status
set -e

# Start the timer
START_TIME=$SECONDS

echo "=========================================="
echo "  refdata-flow Pipeline"
echo "=========================================="
echo ""

RAW_CSV="requests_raw.csv"
OUTDIR="data/references"
DRY_RUN=0


# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --outdir) 
            OUTDIR="$2"
            shift 
            ;;
        --version|-v)
            echo "refdata-flow version 0.9beta20260527"
            exit 0
            ;;
        --dry-run)
            DRY_RUN=1
            ;;
        --clean)
            echo "Cleaning up Nextflow cache, logs, and work directories..."
            rm -rf work/ .nextflow/ .nextflow.log*
            echo "Clean complete! Storage optimized."
            exit 0
            ;;
        *.csv) 
            RAW_CSV="$1" 
            ;;
        *) 
            echo "Unknown parameter passed: $1"
            echo "Usage: bash run_pipeline.sh [requests.csv] [--outdir /path/to/save] [--clean] [--dry-run] [--version]"
            exit 1 

            ;;
    esac
    shift
done

# Ensure genomepy and refgenie are available by creating a local, isolated virtual environment
VENV_DIR=".venv"
if [ ! -d "$VENV_DIR" ]; then
    echo "[Setup] Creating an isolated Python environment for the resolver..."
    python3 -m venv "$VENV_DIR"
    
    echo "[Setup] Installing genomepy and refgenie (this only happens once)..."
    "$VENV_DIR/bin/pip" install --quiet --upgrade pip
    "$VENV_DIR/bin/pip" install --quiet genomepy refgenie
    echo "[Setup] Environment ready!"
    echo ""
fi

if [ ! -f "$RAW_CSV" ]; then
    echo "Error: Cannot find $RAW_CSV!"
    echo "Usage: bash run_pipeline.sh [requests.csv] [--outdir /path/to/save] [--clean] [--dry-run] [--version]"
    exit 1
fi

echo "[Step 1] Running interactive reference resolver..."
"$VENV_DIR/bin/python" bin/resolve.py "$RAW_CSV" requests_resolved.csv

if [ "$DRY_RUN" -eq 1 ]; then
    echo ""
    echo "=========================================="
    echo "  DRY RUN COMPLETE"
    echo "  The following datasets would be downloaded:"
    echo "=========================================="
    "$VENV_DIR/bin/python" -c "import csv; [print(f'  - {r[\"species\"]}: {r[\"assembly\"]} (Provider: {r[\"provider\"]}, Annotations: {r[\"annotation\"]})') for r in csv.DictReader(open('requests_resolved.csv'))]"
    echo "=========================================="
    echo "  Exiting without downloading."
    exit 0
fi

echo ""
echo "[Step 2] Launching Nextflow pipeline..."
# The -resume flag ensures we don't redownload existing datasets
nextflow run main.nf -profile conda --outdir "$OUTDIR" -resume

echo ""
echo "[Step 3] Building Refgenie Configuration..."
"$VENV_DIR/bin/python" bin/update_refgenie.py "$RAW_CSV" requests_resolved.csv "$OUTDIR"

echo ""
ELAPSED_TIME=$(($SECONDS - $START_TIME))
MINUTES=$(($ELAPSED_TIME / 60))
REMAINDER=$(($ELAPSED_TIME % 60))

echo "=========================================="
echo "  Pipeline execution complete!"
echo "  Total runtime: ${MINUTES}m ${REMAINDER}s"
echo "  You can find your new refgenie config at: $OUTDIR/refgenie_export.yaml"
echo "=========================================="
