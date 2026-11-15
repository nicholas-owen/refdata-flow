#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
File a fetched database in the tree beside the vault, and rewrite its catalogue.

    <outdir>/databases/<name>/<version>/     the file as published, or the
                                             contents of the archive it came as
    <outdir>/databases/<name>/<version>/provenance.json
    <outdir>/databases/databases.tsv         one line per database and version

refgenie cannot hold these: every asset there hangs off a genome's digest, and
Swiss-Prot has no genome. So they get the plainest structure that still answers
"which version is this, where did it come from, and has it changed": a directory
per version, a record beside the data, and a table of the lot.

Three rules, the same ones the vault's own index follows:

  * A version is put in place in one step. It is assembled in a hidden directory
    beside its destination and renamed, so a reader never finds half a database,
    and a run that is interrupted leaves nothing that looks finished.
  * A version that is already there is left alone. Its record says what it is;
    a second copy of the same release is not an update.
  * The catalogue is rewritten from the tree on every run, never appended to, so
    it cannot drift from what is on disk when something is removed by hand.

An archive (`unpack: tar` in the definition) is unpacked here, and members that
would land outside the version's directory refuse the whole database. The
checksums in the record are of the file as published, packed.
"""
import argparse
import csv
import json
import os
import shutil
import sys
import tarfile
import tempfile

CATALOGUE_FIELDS = ["database", "version", "path", "file", "file_bytes", "file_sha256",
                    "unpacked", "url", "terms", "retrieved_utc", "checksum_verified"]


def tree(outdir):
    return os.path.join(outdir, "databases")


def held_records(outdir):
    """[(name, version, record)] for every version in the tree that has a record."""
    found = []
    base = tree(outdir)
    if not os.path.isdir(base):
        return found
    for name in sorted(os.listdir(base)):
        if name.startswith(".") or not os.path.isdir(os.path.join(base, name)):
            continue
        for version in sorted(os.listdir(os.path.join(base, name))):
            path = os.path.join(base, name, version, "provenance.json")
            if version.startswith("."):
                continue
            try:
                with open(path) as fh:
                    found.append((name, version, json.load(fh)))
            except (OSError, ValueError):
                continue
    return found


def write_catalogue(outdir):
    """Rewrite <outdir>/databases/databases.tsv from the tree. Returns the line count."""
    base = tree(outdir)
    os.makedirs(base, exist_ok=True)
    rows = [{
        "database": name,
        "version": version,
        "path": f"databases/{name}/{version}/",
        "file": record.get("file", ""),
        "file_bytes": record.get("file_bytes", ""),
        "file_sha256": record.get("file_sha256", ""),
        "unpacked": "yes" if record.get("unpacked") else "no",
        "url": record.get("url", ""),
        "terms": record.get("terms", ""),
        "retrieved_utc": record.get("retrieved_utc", ""),
        "checksum_verified": record.get("checksum_verified", ""),
    } for name, version, record in held_records(outdir)]
    fd, tmp_path = tempfile.mkstemp(dir=base, prefix=".databases-", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CATALOGUE_FIELDS, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(tmp_path, os.path.join(base, "databases.tsv"))
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    return len(rows)


def unpack_tar(archive_path, destination):
    """
    Unpack an archive into a directory, refusing one that reaches outside it.
    Returns the number of files written.
    """
    with tarfile.open(archive_path, "r:*") as archive:
        members = archive.getmembers()
        for member in members:
            parts = [part for part in member.name.split("/") if part not in ("", ".")]
            if member.name.startswith("/") or ".." in parts:
                raise ValueError(f"member '{member.name}' would be written outside the "
                                 "database's directory")
            if member.issym() or member.islnk():
                # A symbolic link's target is relative to the link; a hard link's
                # is a member name, relative to the archive's root.
                base = os.path.dirname(member.name) if member.issym() else ""
                target = os.path.normpath(os.path.join(base, member.linkname))
                if member.linkname.startswith("/") or target.startswith(".."):
                    raise ValueError(f"link '{member.name}' points outside the database's "
                                     "directory")
            elif not (member.isfile() or member.isdir()):
                raise ValueError(f"member '{member.name}' is neither a file, a directory "
                                 "nor a link")
        # 'data' is the standard library's own filter for exactly this (Python 3.12
        # and later); the checks above are what protect an older interpreter.
        if hasattr(tarfile, "data_filter"):
            archive.extractall(destination, filter="data")
        else:
            archive.extractall(destination)
        return sum(1 for member in members if member.isfile())


def main():
    parser = argparse.ArgumentParser(
        description="File a fetched database under <outdir>/databases/ and rewrite the catalogue."
    )
    parser.add_argument("outdir", help="Output directory holding the vault")
    parser.add_argument("--name", help="The database, as its definition names it")
    parser.add_argument("--fetched", help="Directory written by fetch_database.py")
    parser.add_argument("--catalogue-only", action="store_true",
                        help="Only rewrite databases.tsv (and the RO-Crate) from the tree")
    args = parser.parse_args()

    if not args.catalogue_only:
        if not args.name or not args.fetched:
            parser.error("--name and --fetched are required unless --catalogue-only is given")
        with open(os.path.join(args.fetched, f"{args.name}.fetch.json")) as fh:
            record = json.load(fh)
        version = str(record["version"])
        source = os.path.join(args.fetched, record["file"])
        parent = os.path.join(tree(args.outdir), args.name)
        destination = os.path.join(parent, version)
        label = f"{args.name}:{version}"

        if os.path.isfile(os.path.join(destination, "provenance.json")):
            print(f"Database {label} is already in the tree; left as it is.", flush=True)
        else:
            os.makedirs(parent, exist_ok=True)
            # A directory without a record is what an interrupted run leaves.
            shutil.rmtree(destination, ignore_errors=True)
            staging = tempfile.mkdtemp(dir=parent, prefix=f".{version}-")
            try:
                if record.get("unpack") == "tar":
                    try:
                        record["files_unpacked"] = unpack_tar(source, staging)
                    except (tarfile.TarError, EOFError, ValueError) as e:
                        print(f"ERROR: {label} was NOT filed: {e}. Nothing was added.",
                              file=sys.stderr, flush=True)
                        sys.exit(1)
                    record["unpacked"] = True
                else:
                    # A copy, never a link: the tree is what is kept, and it must
                    # not share its bytes with a work directory that is not.
                    shutil.copyfile(source, os.path.join(staging, record["file"]))
                    record["unpacked"] = False
                with open(os.path.join(staging, "provenance.json"), "w") as fh:
                    json.dump(record, fh, indent=4)
                # mkdtemp makes a directory only its owner can enter; give it what
                # any other directory made by this run would have.
                umask = os.umask(0)
                os.umask(umask)
                os.chmod(staging, 0o777 & ~umask)
                os.rename(staging, destination)
            finally:
                shutil.rmtree(staging, ignore_errors=True)
            print(f"Database {label} filed at {destination}"
                  + (f" ({record['files_unpacked']:,} files unpacked)."
                     if record.get("unpacked") else "."), flush=True)

    count = write_catalogue(args.outdir)
    print(f"Wrote the database catalogue ({count} line(s)): "
          f"{os.path.join(tree(args.outdir), 'databases.tsv')}", flush=True)

    # The crate describes everything under the output directory, so it is out of
    # date the moment the tree holds one more database. It needs the vault's own
    # modules and PyYAML, which this step otherwise does not.
    try:
        import update_refgenie as vault
        import vault_crate

        logger = vault.setup_logger(args.outdir)
        vault_crate.write_vault_crate(
            args.outdir, os.path.join(args.outdir, "refgenie_export.yaml"), logger)
    except Exception as e:
        print(f"ERROR: the database is filed, but the vault's RO-Crate could not be "
              f"rewritten: {e}", file=sys.stderr, flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
