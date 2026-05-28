#!/usr/bin/env python3
import argparse
import sys
import json
import os
from datetime import datetime
import genomepy

def main():
    parser = argparse.ArgumentParser(description="Wrapper for genomepy for Nextflow pipeline.")
    parser.add_argument("--species", required=True, help="Species name (e.g., Homo_sapiens)")
    parser.add_argument("--assembly", required=True, help="Assembly name (e.g., GRCh38)")
    parser.add_argument("--provider", required=True, help="Provider (e.g., ensembl, ncbi, ucsc)")
    parser.add_argument("--annotation", action="store_true", help="Download annotation")
    
    args = parser.parse_args()

    # We download into the current working directory.
    # genomepy will automatically create a directory named after the assembly (e.g. GRCh38).
    local_dir = "."
    
    print(f"Downloading {args.assembly} from {args.provider}...")
    
    kwargs = {
        "name": args.assembly,
        "provider": args.provider,
        "genomes_dir": local_dir,
        "annotation": args.annotation,
        "bgzip": True
    }

    try:
        # Install the genome using genomepy's Python API
        genomepy.install_genome(**kwargs)
        
        # Provenance tracking: capture metadata of the download session
        provenance = {
            "timestamp": datetime.now().isoformat(),
            "species": args.species,
            "assembly": args.assembly,
            "provider": args.provider,
            "annotation_requested": args.annotation,
            "tool": "genomepy"
        }
        
        # Write the provenance.json into the created assembly directory
        prov_path = os.path.join(args.assembly, "provenance.json")
        with open(prov_path, "w") as f:
            json.dump(provenance, f, indent=4)
            
        print("Download complete.")
        
    except Exception as e:
        print(f"Error downloading genome: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
