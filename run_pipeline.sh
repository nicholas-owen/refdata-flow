#!/usr/bin/env bash
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later

# Exit immediately if a command exits with a non-zero status
set -e

# Start the timer
START_TIME=$SECONDS

echo "=========================================="
echo "  refdata-flow Pipeline"
echo "=========================================="
echo ""

VERSION="0.9.1"

# Request files live in requests/ so the repo root stays readable, and so the two
# halves of a request - what you asked for, and which assemblies that resolved to -
# sit together. The resolved CSV is not a throwaway build product: bin/resolve.py
# reads it back as its cache, and it records the one scientific choice in the
# pipeline (which assembly), so it must not live anywhere a cleanup could sweep it.
#
# Neither is tracked. RAW_CSV is yours to edit and would otherwise collide with
# every `git pull`; TEMPLATE_CSV is the tracked starting point copied into place on
# a first run. See .gitignore.
RAW_CSV="requests/requests_raw.csv"
RESOLVED_CSV="requests/requests_resolved.csv"
TEMPLATE_CSV="requests/requests_raw.csv.example"
OUTDIR="data/references"
DRY_RUN=0
DEBUG=0
CLEANUP=0
NONINTERACTIVE=0

USAGE="Usage: bash run_pipeline.sh [requests.csv] [--outdir /path/to/save] [--cleanup] [--non-interactive] [--clean] [--dry-run] [--debug] [--version]"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --outdir)
            OUTDIR="$2"
            shift
            ;;
        --version|-v)
            echo "refdata-flow version $VERSION"
            exit 0
            ;;
        --dry-run)
            DRY_RUN=1
            ;;
        --debug)
            DEBUG=1
            ;;
        --non-interactive)
            # Never prompt during resolution. A query matching more than one assembly
            # is reported as unresolved instead, so an unattended run fails rather
            # than hanging on stdin. The pipeline never auto-picks a genome.
            NONINTERACTIVE=1
            ;;
        --cleanup)
            # Remove each downloaded source directory once its contents are verified
            # into the refgenie vault. Opt-in: upstream providers archive old releases,
            # so a deleted download is not guaranteed to be re-fetchable.
            CLEANUP=1
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
            echo "$USAGE"
            exit 1

            ;;
    esac
    shift
done

# refdata-flow pins its Python toolchain in requirements.txt (hash-locked) for
# reproducibility. refgenie 0.13.0 requires Python >= 3.10, so fail fast with a
# clear message rather than a confusing pip resolver error.
REQUIREMENTS="requirements.txt"
if ! python3 -c 'import sys; sys.exit(0 if sys.version_info[:2] >= (3, 10) else 1)'; then
    echo "Error: Python >= 3.10 is required (refgenie 0.13.0)."
    echo "       Found: $(python3 --version 2>&1). Load/activate a newer Python and retry."
    exit 1
fi
if [ ! -f "$REQUIREMENTS" ]; then
    echo "Error: Cannot find $REQUIREMENTS in $(pwd)."
    echo "       This file pins the exact, hash-locked dependency versions."
    exit 1
fi

# Ensure the pinned toolchain is available in a local, isolated virtual environment.
# A marker file records a *successful* install so an interrupted setup is retried
# instead of being silently reused in a broken state.
VENV_DIR=".venv"
VENV_READY="$VENV_DIR/.refdata-flow-install-ok"
if [ ! -f "$VENV_READY" ]; then
    echo "[Setup] Creating an isolated Python environment from $REQUIREMENTS..."
    [ -d "$VENV_DIR" ] || python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --quiet --upgrade pip
    # --require-hashes installs the exact, verified artifacts from the lockfile.
    # genomepy = download engine, refgenie = vault/hashing, pysam = bundled htslib
    # behind the samtools/bgzip shims below, pyyaml = used by update_refgenie.py.
    "$VENV_DIR/bin/pip" install --quiet --require-hashes -r "$REQUIREMENTS"
    touch "$VENV_READY"
    echo "[Setup] Environment ready!"
    echo ""
fi

# ---------------------------------------------------------------------------
# Provide `samtools` and `bgzip` executables inside the venv.
#
# refgenie's fasta recipe shells out to `samtools faidx`, and genomepy shells
# out to `bgzip <file>`. Modern pysam wheels bundle htslib but no longer install
# these CLIs, so we generate thin shims in .venv/bin that dispatch to pysam's
# bundled htslib. The shims call the python interpreter sitting next to them, so
# they keep working even if the project directory is moved. Regenerated on every
# run so a pre-existing .venv from an older version is upgraded transparently.
# ---------------------------------------------------------------------------
"$VENV_DIR/bin/pip" show pysam >/dev/null 2>&1 || "$VENV_DIR/bin/pip" install --quiet --require-hashes -r "$REQUIREMENTS"

cat > "$VENV_DIR/bin/_refdataflow_htsshim.py" <<'PYSHIM_EOF'
#!/usr/bin/env python3
"""samtools / bgzip CLI shims backed by pysam's bundled htslib.

Generated by run_pipeline.sh so a pip-only virtualenv provides the `samtools`
and `bgzip` executables that `refgenie build` (samtools faidx) and genomepy
(bgzip <file>) expect on PATH. Only the command forms used by this pipeline
are implemented."""
import io
import os
import sys


def die(msg, code=1):
    sys.stderr.write("htsshim: " + msg.rstrip() + "\n")
    sys.exit(code)


def main():
    if len(sys.argv) < 2:
        die("missing tool name")
    tool, args = sys.argv[1], sys.argv[2:]
    try:
        import pysam  # noqa: F401
    except ImportError:
        die("pysam is not installed; cannot provide '%s'" % tool)
    if tool == "samtools":
        run_samtools(args)
    elif tool == "bgzip":
        run_bgzip(args)
    else:
        die("unknown tool '%s'" % tool)


def run_samtools(args):
    import pysam
    if not args:
        die("samtools: no subcommand given")
    sub, rest = args[0], args[1:]
    disp = getattr(getattr(pysam, "samtools", None), sub, None)
    if disp is None:
        die("samtools: subcommand '%s' not supported by shim" % sub)
    try:
        # catch_stdout=False -> stream to the real stdout, like the CLI.
        disp(*rest, catch_stdout=False)
    except Exception as e:  # pysam raises SamtoolsError on non-zero exit
        die("samtools %s failed: %s" % (sub, e))


def run_bgzip(argv):
    import pysam
    decompress = stdout = force = keep = False
    files = []
    it = iter(argv)
    for a in it:
        if a in ("-d", "--decompress"):
            decompress = True
        elif a in ("-c", "--stdout"):
            stdout = True
        elif a in ("-f", "--force"):
            force = True
        elif a in ("-k", "--keep"):
            keep = True
        elif a in ("-@", "--threads"):
            next(it, None)            # consume + ignore the thread count
        elif a.startswith("-@"):
            pass                       # -@4 form
        elif a.startswith("-"):
            pass                       # ignore other flags rather than crash
        else:
            files.append(a)

    if not files:                      # stream mode: stdin -> stdout
        data = sys.stdin.buffer.read()
        if decompress:
            bf = pysam.BGZFile(fileobj=io.BytesIO(data), mode="rb")
            sys.stdout.buffer.write(bf.read()); bf.close()
        else:
            buf = io.BytesIO(); bf = pysam.BGZFile(fileobj=buf, mode="wb")
            bf.write(data); bf.close(); sys.stdout.buffer.write(buf.getvalue())
        return

    for f in files:
        if decompress:
            out = f[:-3] if f.endswith(".gz") else f + ".decompressed"
            bf = pysam.BGZFile(f, "rb"); data = bf.read(); bf.close()
            if stdout:
                sys.stdout.buffer.write(data)
            else:
                if os.path.exists(out) and not force:
                    die("%s already exists; use -f to overwrite" % out)
                with open(out, "wb") as fh:
                    fh.write(data)
                if not keep:
                    os.remove(f)
        else:
            out = f + ".gz"
            if stdout:
                buf = io.BytesIO(); bf = pysam.BGZFile(fileobj=buf, mode="wb")
                with open(f, "rb") as fh:
                    bf.write(fh.read())
                bf.close(); sys.stdout.buffer.write(buf.getvalue())
            else:
                if os.path.exists(out) and not force:
                    die("%s already exists; use -f to overwrite" % out)
                pysam.tabix_compress(f, out, force=force)
                if not keep:
                    os.remove(f)


if __name__ == "__main__":
    main()
PYSHIM_EOF

cat > "$VENV_DIR/bin/samtools" <<'SAMTOOLS_EOF'
#!/bin/sh
# pysam-backed samtools shim (refdata-flow). Uses the python next to this script.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$DIR/python" "$DIR/_refdataflow_htsshim.py" samtools "$@"
SAMTOOLS_EOF

cat > "$VENV_DIR/bin/bgzip" <<'BGZIP_EOF'
#!/bin/sh
# pysam-backed bgzip shim (refdata-flow). Uses the python next to this script.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$DIR/python" "$DIR/_refdataflow_htsshim.py" bgzip "$@"
BGZIP_EOF

chmod +x "$VENV_DIR/bin/samtools" "$VENV_DIR/bin/bgzip"

# Put the venv bin on PATH so refgenie's `samtools faidx` subprocess (and any
# genomepy `bgzip` call) resolves to the shims above.
VENV_BIN="$(cd "$VENV_DIR/bin" && pwd)"
export PATH="$VENV_BIN:$PATH"

if [ ! -f "$RAW_CSV" ]; then
    # First run: seed the request file from the tracked template, then stop. Copying
    # and carrying straight on would download whatever the template happens to list,
    # which is not what someone who has not written a request yet is asking for.
    if [ "$RAW_CSV" = "requests/requests_raw.csv" ] && [ -f "$TEMPLATE_CSV" ]; then
        mkdir -p requests
        cp "$TEMPLATE_CSV" "$RAW_CSV"
        echo "=========================================="
        echo "  Created $RAW_CSV from the template."
        echo ""
        echo "  Edit it to list the genomes you want, then re-run. It is not"
        echo "  tracked by git, so your list will not collide with updates."
        echo "=========================================="
        exit 0
    fi
    echo "Error: Cannot find $RAW_CSV!"
    # Requests moved into requests/ after 0.9.1. Say so rather than reporting a bare
    # not-found, and do not silently read the old path - that would make the run do
    # something other than what the docs describe.
    if [ -f "requests_raw.csv" ]; then
        echo "Note: found requests_raw.csv in the repo root. Request files now live"
        echo "      in requests/. Move it with: mkdir -p requests && git mv requests_raw.csv $RAW_CSV"
    fi
    echo "$USAGE"
    exit 1
fi

echo "[Step 1] Running interactive reference resolver..."
# resolve.py exits 2 when some rows could not be resolved (as opposed to 1 for a
# fatal error). Those rows used to be skipped silently while the run continued, so
# a user could ask for five genomes, receive three, and never be told. Stop here
# instead: nothing has been downloaded yet, and the request CSV is the thing to fix.
set +e
RESOLVE_ARGS=("$RAW_CSV" "$RESOLVED_CSV")
if [ "$NONINTERACTIVE" -eq 1 ]; then
    RESOLVE_ARGS+=(--non-interactive)
fi
"$VENV_DIR/bin/python" bin/resolve.py "${RESOLVE_ARGS[@]}"
RESOLVE_RC=$?
set -e
if [ "$RESOLVE_RC" -eq 2 ]; then
    echo ""
    echo "=========================================="
    echo "  Stopping: some requests could not be resolved (see above)."
    echo "  Fix or remove those rows in $RAW_CSV and re-run."
    echo "  Requests that did resolve have been saved to $RESOLVED_CSV."
    echo "=========================================="
    exit 2
elif [ "$RESOLVE_RC" -ne 0 ]; then
    exit "$RESOLVE_RC"
fi

if [ "$DRY_RUN" -eq 1 ]; then
    echo ""
    echo "=========================================="
    echo "  DRY RUN COMPLETE"
    echo "  The following datasets would be downloaded:"
    echo "=========================================="
    # Release is shown because pinning it is the point of the dry run: 'current'
    # means whatever the provider publishes on the day the download actually runs,
    # which is not reproducible. .get() rather than [] so a resolved CSV written
    # before the release column existed still prints.
    #
    # The path comes in as sys.argv[1], not interpolated into this string: the
    # snippet is already three levels of nested quoting deep, and a shell variable
    # inside it is how that becomes unreadable and then wrong.
    "$VENV_DIR/bin/python" -c "import csv, sys; [print(f'  - {r[\"species\"]}: {r[\"assembly\"]} (Provider: {r[\"provider\"]}, Annotations: {r[\"annotation\"]}, Release: {r.get(\"release\") or \"current\"})') for r in csv.DictReader(open(sys.argv[1]))]" "$RESOLVED_CSV"
    echo "=========================================="
    echo "  Exiting without downloading."
    exit 0
fi

echo ""
echo "[Step 2] Launching Nextflow pipeline..."
# Assemble Nextflow arguments. --debug streams each task's progress messages to
# the console (process `debug`); -ansi-log false makes that streamed output
# readable instead of being repainted over by the single-line status display.
# --input is passed explicitly rather than left to nextflow.config's default. The
# wrapper and the config each held their own copy of this path, agreeing only by
# coincidence; passing it here makes the wrapper authoritative and leaves the config
# default for a bare `nextflow run main.nf`.
NF_ARGS="-profile conda --input $RESOLVED_CSV --outdir $OUTDIR -resume"
if [ "$DEBUG" -eq 1 ]; then
    echo "[Step 2] Debug mode: streaming per-task download progress."
    NF_ARGS="$NF_ARGS --debug true -ansi-log false"
fi
# The -resume flag ensures we don't redownload existing datasets
nextflow run main.nf $NF_ARGS

echo ""
echo "[Step 3] Building Refgenie Configuration..."
REFGENIE_ARGS=("$RAW_CSV" "$RESOLVED_CSV" "$OUTDIR")
if [ "$CLEANUP" -eq 1 ]; then
    echo "[Step 3] --cleanup: verified source directories will be removed after ingest."
    REFGENIE_ARGS+=(--cleanup)
fi
"$VENV_DIR/bin/python" bin/update_refgenie.py "${REFGENIE_ARGS[@]}"

echo ""
ELAPSED_TIME=$(($SECONDS - $START_TIME))
MINUTES=$(($ELAPSED_TIME / 60))
REMAINDER=$(($ELAPSED_TIME % 60))

echo "=========================================="
echo "  Pipeline execution complete!"
echo "  Total runtime: ${MINUTES}m ${REMAINDER}s"
echo "  You can find your new refgenie config at: $OUTDIR/refgenie_export.yaml"
if [ "$CLEANUP" -eq 0 ]; then
    echo ""
    echo "  Downloaded sources were retained under $OUTDIR/<provider>/<species>/."
    echo "  Re-run with --cleanup to remove them once you are satisfied with the vault."
fi
echo "=========================================="
# end of run_pipeline.sh
