#!/usr/bin/env python3
# Copyright (C) 2026 Nicholas Owen
# SPDX-License-Identifier: GPL-3.0-or-later
"""
Register one derived asset in the refgenie vault, against a genome already there.

bin/update_refgenie.py catalogues what was DOWNLOADED: a genome and its
annotation. This script catalogues what was BUILT from them afterwards by a derive
process - a transcriptome today. The split mirrors the workflow: deriving is
heavy and runs wherever the profile sends it; registering is a few file copies
that must happen one at a time, on the machine that holds the vault.

It does five things, and the last four are the reason it is not a bare
`refgenie build`:

  1. Builds <genome>/<asset>:<tag> from the files it is given, with a recipe from
     recipes/, naming the parent assets so refgenie records the lineage.
  2. Confirms the asset's files are actually in the vault. A skipped build still
     writes its seek keys to the config (see update_refgenie.fasta_has_dictionary),
     so the config alone cannot be believed.
  3. Moves the default tag only forwards, by the same rule as the annotation it
     was derived from (update_refgenie.should_become_default).
  4. Refreshes the alias links, so the asset is reachable under every name the
     genome has and not only the one it was built under.
  5. Rewrites the vault's index and its RO-Crate, which describe everything the
     vault holds and so are out of date the moment it holds one more thing.
"""
import argparse
import os
import subprocess
import sys

import yaml

import update_refgenie as vault
import vault_crate


def find_digest(config, genome):
    for digest, data in (config.get("genomes") or {}).items():
        if genome in ((data or {}).get("aliases") or []):
            return digest
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Register a derived asset against a genome already in the vault."
    )
    parser.add_argument("outdir", help="Output directory holding the refgenie vault")
    parser.add_argument("--genome", required=True,
                        help="The genome's vault name, e.g. GRCh38_p14_ensembl")
    parser.add_argument("--asset", required=True,
                        help="Asset to build; recipes/<asset>.json is its recipe")
    parser.add_argument("--tag", required=True, help="Tag to build the asset under")
    parser.add_argument("--file", action="append", default=[], metavar="KEY=PATH",
                        help="An input file the recipe requires. Repeatable.")
    parser.add_argument("--parent", action="append", default=[], metavar="KEY=ASSET[:TAG]",
                        help="A parent asset of this genome, as the recipe names it, "
                             "e.g. annotation=ensembl_gtf:116. Repeatable.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    logger = vault.setup_logger(args.outdir)
    config_file = os.path.join(args.outdir, "refgenie_export.yaml")
    refgenie_bin = os.path.join(os.path.dirname(sys.executable), "refgenie")
    target = f"{args.genome}/{args.asset}:{args.tag}"

    recipe = os.path.join(vault.RECIPES_DIR, f"{args.asset}.json")
    if not os.path.isfile(recipe):
        logger.error("No recipe for asset '%s' (looked for '%s')", args.asset,
                     os.path.abspath(recipe))
        sys.exit(1)

    try:
        with open(config_file, "r") as f:
            config = yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError) as e:
        logger.error("Cannot read the vault config '%s': %s", config_file, e)
        sys.exit(1)
    if find_digest(config, args.genome) is None:
        logger.error("Genome '%s' is not in the vault, so '%s' cannot be registered "
                     "against it", args.genome, target)
        sys.exit(1)

    files = []
    for item in args.file:
        key, separator, path = item.partition("=")
        if not separator or not key or not os.path.isfile(path):
            logger.error("--file expects KEY=PATH to an existing file, got '%s'", item)
            sys.exit(1)
        files.append(f"{key}={os.path.abspath(path)}")
    parents = []
    for item in args.parent:
        key, separator, asset = item.partition("=")
        if not separator or not key or not asset:
            logger.error("--parent expects KEY=ASSET[:TAG], got '%s'", item)
            sys.exit(1)
        parents.append(f"{key}={args.genome}/{asset}")

    # Read BEFORE the build, for the same reason as in update_refgenie: afterwards
    # there is no telling whether refgenie made the new tag the default itself.
    previous_default = vault.current_default_tag(config_file, args.genome, args.asset)

    command = [sys.executable, refgenie_bin, "build", target,
               "--recipe", os.path.abspath(recipe), "-c", config_file]
    if files:
        command += ["--files"] + files
    if parents:
        command += ["--assets"] + parents
    logger.info("Building '%s' (recipe: %s [custom]; parents: %s)", target, args.asset,
                ", ".join(parents) or "none given")
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logger.error("Build of '%s' failed: %s", target, e)
        sys.exit(1)

    with open(config_file, "r") as f:
        config = yaml.safe_load(f) or {}
    digest = find_digest(config, args.genome)
    if not vault.asset_file(config, config_file, digest, args.asset, args.tag, args.asset):
        logger.error("'%s' is listed in the config but its file is not in the vault",
                     target)
        sys.exit(1)

    if vault.should_become_default(previous_default, args.tag):
        try:
            subprocess.run([sys.executable, refgenie_bin, "tag", target,
                            "--default", "-c", config_file], check=True)
        except subprocess.CalledProcessError as e:
            logger.warning("Built '%s' but could not make it the default tag: %s. It is "
                           "addressable by tag.", target, e)
    else:
        logger.info("Built '%s'. The default tag stays at '%s', which is a newer "
                    "release; name the tag to use this one.", target, previous_default)

    vault.relink_aliases(config_file, args.genome, logger)
    logger.info("Registered '%s'", target)

    # The index and the crate describe the whole vault, so both are rewritten now
    # that it holds one more asset. Left to update_refgenie alone they would be
    # stale for everything registered after the genomes in the same run.
    try:
        vault.write_vault_index(args.outdir, config_file, logger)
        vault_crate.write_vault_crate(args.outdir, config_file, logger)
    except Exception as e:
        logger.error("Registered '%s', but could not rewrite the vault's index and "
                     "RO-Crate: %s", target, e)
        sys.exit(1)


if __name__ == "__main__":
    main()
