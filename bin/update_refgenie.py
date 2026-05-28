#!/usr/bin/env python3
import sys
import os
import csv
import glob
import subprocess
import shutil

def main():
    if len(sys.argv) < 4:
        print("Usage: python bin/update_refgenie.py <requests_raw.csv> <requests_resolved.csv> <outdir>")
        sys.exit(1)

    raw_csv = sys.argv[1]
    resolved_csv = sys.argv[2]
    outdir = sys.argv[3]
    
    config_file = os.path.join(outdir, "refgenie_export.yaml")
    refgenie_bin = os.path.join(os.path.dirname(sys.executable), "refgenie")

    # Ensure refgenie config is initialized
    if not os.path.exists(config_file):
        print(f"Initializing refgenie config at {config_file}...")
        subprocess.run([sys.executable, refgenie_bin, "init", "-c", config_file], check=True)

    # Map generic raw queries to exact species
    query_map = {}
    custom_aliases_map = {}
    try:
        with open(raw_csv, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                sp = row.get("species", "").strip()
                q = row.get("query", "").strip()
                aliases_str = row.get("aliases", "").strip()
                
                if sp:
                    if q:
                        query_map[sp] = q
                    if aliases_str:
                        # Split by semicolon and strip whitespace to allow multiple aliases
                        custom_aliases_map[sp] = [a.strip() for a in aliases_str.split(";") if a.strip()]
    except Exception as e:
        print(f"Warning: Could not read {raw_csv}: {e}")

    # Process each resolved request
    try:
        with open(resolved_csv, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                species = row.get("species", "").strip()
                assembly = row.get("assembly", "").strip()
                provider = row.get("provider", "").strip()
                annotation = row.get("annotation", "false").strip().lower() == "true"
                
                print(f"\nProcessing {assembly}...")
                
                asm_dir = os.path.join(outdir, provider, species, assembly)
                assembly_safe = assembly.replace(".", "_")

                if os.path.exists(asm_dir):
                    # Find fasta
                    fasta_files = glob.glob(os.path.join(asm_dir, "*.fa")) + glob.glob(os.path.join(asm_dir, "*.fa.gz"))
                    if fasta_files:
                        fasta_path = os.path.abspath(fasta_files[0])
                        print(f"  Building FASTA digest in refgenie for {assembly_safe}...")
                        cmd = [sys.executable, refgenie_bin, "build", f"{assembly_safe}/fasta", "--files", f"fasta={fasta_path}", "-c", config_file]
                        subprocess.run(cmd, check=True)
                    else:
                        print(f"  Warning: No FASTA found for {assembly}")

                    # Find annotation
                    if annotation:
                        gtf_files = glob.glob(os.path.join(asm_dir, "*.gtf")) + glob.glob(os.path.join(asm_dir, "*.gtf.gz"))
                        if gtf_files:
                            gtf_path = os.path.abspath(gtf_files[0])
                            print(f"  Building GTF annotation digest in refgenie for {assembly_safe}...")
                            
                            recipe_name = "ensembl_gtf" if provider.lower() == "ensembl" else "gencode_gtf"
                            cmd = [sys.executable, refgenie_bin, "build", f"{assembly_safe}/{recipe_name}", "--files", f"{recipe_name}={gtf_path}", "-c", config_file]
                            subprocess.run(cmd, check=False)

                    # Clean Up the entire Nextflow Provider directory tree safely
                    print(f"  Cleaning up Nextflow temporary download directory to save space: {asm_dir}")
                    shutil.rmtree(asm_dir, ignore_errors=True)
                else:
                    print(f"  Directory {asm_dir} not found. Assuming already in Refgenie vault.")

                # -------------------------------------------------------------
                # ALIAS ASSIGNMENT (Direct YAML update to avoid CLI/OS prompts)
                # -------------------------------------------------------------
                aliases_to_set = set([assembly_safe, assembly, assembly.split(".")[0]])
                
                raw_query = query_map.get(species)
                if raw_query:
                    aliases_to_set.add(raw_query)
                    
                for ca in custom_aliases_map.get(species, []):
                    aliases_to_set.add(ca)

                try:
                    import yaml
                    with open(config_file, "r") as f:
                        config = yaml.safe_load(f)
                    
                    target_digest = None
                    if "genomes" in config:
                        for digest, data in config["genomes"].items():
                            existing_aliases = data.get("aliases", [])
                            # If ANY of our desired aliases are currently pointing to this digest
                            if any(a in existing_aliases for a in aliases_to_set):
                                target_digest = digest
                                break
                    
                    if target_digest:
                        # Merge aliases and sort for clean diffs
                        current_aliases = config["genomes"][target_digest].get("aliases", [])
                        merged_aliases = sorted(list(set(current_aliases) | aliases_to_set))
                        config["genomes"][target_digest]["aliases"] = merged_aliases
                        
                        print(f"  Applying {len(merged_aliases)} aliases to digest {target_digest[:8]}...: {merged_aliases}")
                        with open(config_file, "w") as f:
                            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
                    else:
                        print(f"  Warning: Could not resolve a digest for {assembly} in config. Skipping alias update.")
                except Exception as e:
                    print(f"  Failed to update aliases via YAML: {e}")
                
    except Exception as e:
        print(f"Error processing resolved requests: {e}")
        sys.exit(1)
        
    print("\nRefgenie config updated successfully!")

if __name__ == "__main__":
    main()
