#!/usr/bin/env python3
import sys
import os
import csv
import genomepy

def main():
    if len(sys.argv) < 2:
        print("Usage: python bin/resolve.py <requests_raw.csv> [output_resolved.csv]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else "requests_resolved.csv"

    # Load existing resolutions to act as a cache
    resolved_cache = {}
    if os.path.exists(output_file):
        try:
            with open(output_file, mode='r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sp = row.get("species", "").strip()
                    if sp:
                        resolved_cache[sp] = row
        except Exception:
            pass

    resolved_rows = []

    try:
        with open(input_file, mode='r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                species = row.get("species", "").strip()
                query = row.get("query", "").strip()
                provider_filter = row.get("provider", "").strip().lower()
                annotation = row.get("annotation", "false").strip()

                if not query:
                    continue

                # Check if this species was already resolved in a previous run
                if species in resolved_cache:
                    cached_row = resolved_cache[species]
                    print(f"\n--- Species '{species}' already resolved to {cached_row['assembly']}. Skipping search. ---")
                    # Update annotation flag in case they changed it in raw file
                    cached_row['annotation'] = annotation
                    resolved_rows.append(cached_row)
                    continue

                print(f"\n--- Resolving query: '{query}' for {species} ---")
                
                # genomepy search returns a generator of lists
                results = list(genomepy.search(query))
                
                # Filter by provider if the user specified one
                if provider_filter:
                    results = [r for r in results if r[1].lower() == provider_filter]
                
                if not results:
                    msg = f"ERROR: No matches found for '{query}'"
                    if provider_filter:
                        msg += f" with provider '{provider_filter}'"
                    print(msg + ". Skipping.")
                    continue
                
                if len(results) == 1:
                    match = results[0]
                    assembly = match[0]
                    provider = match[1]
                    print(f"Found exact match: {assembly} from {provider}.")
                else:
                    # Multiple matches, interactive prompt
                    print(f"Found {len(results)} matches for '{query}':")
                    for idx, match in enumerate(results):
                        name = match[0]
                        provider = match[1]
                        extra = " | ".join(str(x) for x in match[2:6] if x)
                        print(f"  [{idx + 1}] {name} (Provider: {provider}) | {extra}")
                    
                    choice = 0
                    while choice < 1 or choice > len(results):
                        try:
                            choice_str = input(f"Select an option [1-{len(results)}]: ")
                            choice = int(choice_str)
                        except ValueError:
                            pass
                    
                    selected_match = results[choice - 1]
                    assembly = selected_match[0]
                    provider = selected_match[1]
                    print(f"Selected: {assembly} from {provider}.")

                resolved_rows.append({
                    "species": species,
                    "assembly": assembly,
                    "provider": provider,
                    "annotation": annotation
                })
                
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    with open(output_file, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["species", "assembly", "provider", "annotation"])
        writer.writeheader()
        for row in resolved_rows:
            writer.writerow(row)
            
    print(f"\nSuccessfully resolved {len(resolved_rows)} requests. Written to {output_file}.")

if __name__ == "__main__":
    main()
