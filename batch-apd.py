#!/usr/bin/env python3
# Orchestrates the combined batch processing pipeline for Amyloid Packing Difference (APD)
# Extracts contacts from structural files and performs all-vs-all pairwise comparisons

import os
import sys
import subprocess
import shutil
import itertools
import re
import csv
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_comparison(pair, contacts_dir, logs_dir, figs_dir, script_path):
    """Executes a single pairwise comparison via subprocess."""
    file1, file2 = pair
    
    # Executes compare.py from within the Contacts directory to prevent pathing issues
    cmd = [sys.executable, str(script_path), file1.name, file2.name]
    result = subprocess.run(cmd, cwd=contacts_dir, capture_output=True, text=True)
    
    base1 = file1.stem
    base2 = file2.stem
    log_name = f"{base1}_vs_{base2}.log"
    svg_name = f"{base1}_vs_{base2}.svg"
    
    # Moves the generated output files to their respective organized directories
    src_log = contacts_dir / log_name
    if src_log.exists():
        shutil.move(str(src_log), str(logs_dir / log_name))
        
    src_svg = contacts_dir / svg_name
    if src_svg.exists():
        shutil.move(str(src_svg), str(figs_dir / svg_name))
            
    return file1, file2, result.returncode, logs_dir / log_name

def main():
    """Executes the full batch APD workflow."""
    # Sets up base directories and paths for inputs and outputs
    base_dir = Path.cwd()
    input_dir = base_dir / "PDBs"
    contacts_dir = base_dir / "Contacts"
    comparisons_dir = base_dir / "Comparisons"
    logs_dir = comparisons_dir / "Logs"
    figs_dir = comparisons_dir / "Figures"
    
    # Defaults to the APD submodule script paths and falls back to local paths if needed
    contacts_script = base_dir / "APD/contacts.py"
    compare_script = base_dir / "APD/compare.py"

    if not contacts_script.exists() or not compare_script.exists():
        contacts_script = base_dir / "contacts.py"
        compare_script = base_dir / "compare.py"
        if not contacts_script.exists() or not compare_script.exists():
            print("Error: Could not find contacts.py or compare.py.")
            return

    # Ensures all required directories exist before processing begins
    input_dir.mkdir(parents=True, exist_ok=True)
    contacts_dir.mkdir(parents=True, exist_ok=True)
    comparisons_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    # Initiates Phase 1 to extract contacts from all structure files
    print("=== Phase 1: Extracting Contacts ===")
    structure_files = []
    for ext in ("*.pdb", "*.cif", "*.mmcif"):
        structure_files.extend(input_dir.glob(ext))

    if not structure_files:
        print(f"No structure files found in '{input_dir}'. Please add PDB/CIF files and run again.")
    else:
        print(f"Found {len(structure_files)} structures to process in '{input_dir}'.\n")
        for struct_file in structure_files:
            print(f"Processing: {struct_file.name}...")
            
            result = subprocess.run(
                [sys.executable, str(contacts_script), str(struct_file)],
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                print(f"  -> [!] Failed to process {struct_file.name}")
                print(f"  -> Error details: {result.stderr.strip()}")
                continue

            base_name = struct_file.stem
            expected_csv = struct_file.parent / f"{base_name}_contacts.csv"
            expected_bild = struct_file.parent / f"{base_name}_contacts.bild"

            # Relocates the extracted contacts and bild files to the centralized contacts directory
            for out_file in (expected_csv, expected_bild):
                if out_file.exists():
                    dest_file = contacts_dir / out_file.name
                    shutil.move(str(out_file), str(dest_file))
                    
            print(f"  -> Successfully extracted contacts to '{contacts_dir.name}/'")

    # Initiates Phase 2 to perform all-vs-all pairwise comparisons
    print("\n=== Phase 2: Pairwise Comparisons ===")
    csv_files = sorted(list(contacts_dir.glob("*_contacts.csv")))
    if len(csv_files) < 2:
        print("Not enough contact CSVs in 'Contacts' to perform pairwise comparisons. Exiting.")
        return

    # Generates all unique combinatorial pairs for the pairwise comparisons
    pairs = list(itertools.combinations(csv_files, 2))
    print(f"Found {len(csv_files)} structure CSVs. Generating {len(pairs)} pairwise comparisons...")

    log_files = []
    # Determines the number of available CPU cores to allocate for parallel processing
    cpu_cores = os.cpu_count() or 4
    workers = max(1, cpu_cores - 2)
    print(f"Using {workers} CPU cores for parallel processing...\n")
    
    # Distributes pairwise comparisons across the process pool
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(run_comparison, pair, contacts_dir, logs_dir, figs_dir, compare_script): pair 
            for pair in pairs
        }
        
        for i, future in enumerate(as_completed(futures), 1):
            f1, f2, returncode, log_path = future.result()
            if returncode == 0 and log_path.exists():
                log_files.append((f1.stem, f2.stem, log_path))
            else:
                print(f"  -> [!] Comparison failed: {f1.name} vs {f2.name}")
            
            if i % 10 == 0 or i == len(pairs):
                print(f"  -> Progress: {i}/{len(pairs)} comparisons completed")

    # Initiates Phase 3 to parse logs and generate the distance matrices
    print("\n=== Phase 3: Distance Matrix Generation ===")
    # Extracts unique structure names by removing the contacts suffix
    structure_names = sorted(list({f.stem.replace("_contacts", "") for f in csv_files}))
    
    # Initializes symmetric matrices with zeros on the diagonal
    matrix_xy = {s1: {s2: 0.0 for s2 in structure_names} for s1 in structure_names}
    matrix_xyz = {s1: {s2: 0.0 for s2 in structure_names} for s1 in structure_names}
    
    # Compiles regular expressions to locate the XY and XYZ APD percentage scores within the text logs
    score_pattern_xy = re.compile(r"Amyloid Packing Difference \(XY\): ([\d\.]+)%")
    score_pattern_xyz = re.compile(r"Amyloid Packing Difference \(XYZ\): ([\d\.]+)%")
    
    # Parses each log file to populate the distance matrices
    for stem1, stem2, log_path in log_files:
        name1 = stem1.replace("_contacts", "")
        name2 = stem2.replace("_contacts", "")
        
        with open(log_path, "r") as f:
            content = f.read()
            
        scores_xy = score_pattern_xy.findall(content)
        scores_xyz = score_pattern_xyz.findall(content)
        
        avg_score_xy = sum(float(s) for s in scores_xy) / len(scores_xy) if scores_xy else 100.0
        avg_score_xyz = sum(float(s) for s in scores_xyz) / len(scores_xyz) if scores_xyz else 100.0
            
        matrix_xy[name1][name2] = avg_score_xy
        matrix_xy[name2][name1] = avg_score_xy
        matrix_xyz[name1][name2] = avg_score_xyz
        matrix_xyz[name2][name1] = avg_score_xyz

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

    print("\n[+] Full batch APD workflow completed successfully!")

if __name__ == "__main__":
    main()