#!/usr/bin/env python3
# Orchestrates the combined batch processing pipeline for Amyloid Packing Difference (APD)
# Extracts contacts from structural files and performs all-vs-all pairwise comparisons

import os
import sys

# Prevents Python from generating __pycache__ directories when importing modules
sys.dont_write_bytecode = True
import shutil
import itertools
import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Dynamically appends the APD script directory to the system path across all OS process spawn methods
base_dir = Path.cwd()
apd_dir = base_dir / "APD"
if apd_dir.exists() and (apd_dir / "compare.py").exists():
    sys.path.insert(0, str(apd_dir))
else:
    sys.path.insert(0, str(base_dir))

def process_pair(pair, logs_dir, figs_dir):
    """Executes a single pairwise comparison natively returning the APD scores"""
    file1, file2 = pair
    base1 = file1.stem.replace("_contacts", "")
    base2 = file2.stem.replace("_contacts", "")
    
    log_path = logs_dir / f"{base1}_vs_{base2}.log"
    svg_path = figs_dir / f"{base1}_vs_{base2}.svg"
    
    try:
        import compare
        avg_xy, avg_xyz = compare.run_comparison(str(file1), str(file2), svg_path=str(svg_path), log_path=str(log_path), quiet=True)
        return base1, base2, True, avg_xy, avg_xyz
    except Exception as e:
        return base1, base2, False, 100.0, 100.0

def main():
    """Executes the full batch APD workflow"""
    # Sets up base directories and paths for inputs and outputs
    input_dir = base_dir / "PDBs"
    contacts_dir = base_dir / "Contacts"
    comparisons_dir = base_dir / "Comparisons"
    logs_dir = comparisons_dir / "Logs"
    figs_dir = comparisons_dir / "Figures"
    
    try:
        import compare
    except ImportError:
        print("  -> [!] Error: Could not import compare.py; ensure it is located in the root or APD directory")
        return

    # Ensures all required directories exist before processing begins
    input_dir.mkdir(parents=True, exist_ok=True)
    contacts_dir.mkdir(parents=True, exist_ok=True)
    comparisons_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # Extracts contacts from all structure files
    print("\nExtracting contacts of all structures in PDBs/...\n")
    structure_files = []
    for ext in ("*.pdb", "*.cif", "*.mmcif"):
        structure_files.extend(input_dir.glob(ext))

    if not structure_files:
        print(f"  -> [!] No structure files found in '{input_dir}'; please add PDB/CIF files and run again")
    else:
        print(f"  -> Found {len(structure_files)} structures to process in '{input_dir}'\n")
        
        import contacts
        for struct_file in structure_files:
            print(f"Processing: {struct_file.name}...")
            
            try:
                base_name = struct_file.stem
                expected_csv = contacts_dir / f"{base_name}_contacts.csv"
                expected_bild = contacts_dir / f"{base_name}_contacts.bild"
                
                # Evaluates intelligent caching to skip extraction if output files already exist and are current
                if expected_csv.exists() and expected_bild.exists():
                    in_mtime = struct_file.stat().st_mtime
                    if expected_csv.stat().st_mtime > in_mtime and expected_bild.stat().st_mtime > in_mtime:
                        print(f"  -> Skipping: Contacts already up-to-date in '{contacts_dir.name}/'")
                        continue

                # Executes the modularized contact extraction natively
                contacts.extract_contacts(str(struct_file), csv_file=str(expected_csv), bild_file=str(expected_bild), quiet=True)
                print(f"  -> Successfully extracted contacts to '{contacts_dir.name}/'")
            except Exception as e:
                print(f"  -> [!] Failed to process {struct_file.name}")
                print(f"  -> Error details: {e}")
                continue

    # Performs all-vs-all pairwise comparisons of contacts
    print("\n\nPerforming pairwise comparisons of all structure contact lists in Contacts/...\n")
    csv_files = sorted(contacts_dir.glob("*_contacts.csv"))
    if len(csv_files) < 2:
        print("  -> [!] Not enough contact CSVs in 'Contacts' to perform pairwise comparisons")
        return

    # Generates all unique combinatorial pairs for the pairwise comparisons
    pairs = list(itertools.combinations(csv_files, 2))
    print(f"  -> Found {len(csv_files)} structure CSVs; generating {len(pairs)} pairwise comparisons")

    log_files = []
    # Determines the number of available CPU cores to allocate for parallel processing
    cpu_cores = os.cpu_count() or 4
    workers = max(1, cpu_cores - 1)
    print(f"  -> Using {workers} CPU cores for parallel processing\n")
    
    # Distributes pairwise comparisons across the process pool
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(process_pair, pair, logs_dir, figs_dir): pair 
            for pair in pairs
        }
        
        # Evaluates the completed futures and captures the returned APD scores directly
        for i, future in enumerate(as_completed(futures), 1):
            base1, base2, success, avg_xy, avg_xyz = future.result()
            if success:
                log_files.append((base1, base2, avg_xy, avg_xyz))
            else:
                print(f"  -> [!] Comparison failed internally: {base1} vs {base2}")
            
            if i % 10 == 0 or i == len(pairs):
                print(f"  -> Progress: {i}/{len(pairs)} comparisons completed")

    # Parses logs and generates the distance matrices
    print("\n\nWriting APD difference matrices to CSV files...\n")
    # Extracts unique structure names by removing the contacts suffix
    structure_names = sorted({f.stem.replace("_contacts", "") for f in csv_files})
    
    # Initializes symmetric matrices with zeros on the diagonal
    matrix_xy = {s1: {s2: 0.0 for s2 in structure_names} for s1 in structure_names}
    matrix_xyz = {s1: {s2: 0.0 for s2 in structure_names} for s1 in structure_names}
    
    # Populates the matrices using the directly returned score data
    for name1, name2, avg_xy, avg_xyz in log_files:
        matrix_xy[name1][name2] = avg_xy
        matrix_xy[name2][name1] = avg_xy
        matrix_xyz[name1][name2] = avg_xyz
        matrix_xyz[name2][name1] = avg_xyz

    # Writes the fully populated distance matrices to CSV files for downstream plotting
    for label, matrix in [("XY", matrix_xy), ("XYZ", matrix_xyz)]:
        matrix_out = comparisons_dir / f"apd_distance_matrix_{label}.csv"
        with open(matrix_out, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Structure"] + structure_names)
            for s1 in structure_names:
                row = [s1] + [f"{matrix[s1][s2]:.1f}" for s2 in structure_names]
                writer.writerow(row)
        print(f"  -> Distance matrix ({label}) written to '{matrix_out}'")

    print("\n\nFull batch APD workflow completed successfully!")

if __name__ == "__main__":
    main()