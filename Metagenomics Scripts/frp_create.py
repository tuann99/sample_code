#!/usr/bin/env python3

import os
import re
import argparse
import matplotlib.pyplot as plt
from simplesam import Reader
from tqdm import tqdm

def calculate_percent_identity(md_field):
    matches = sum(map(int, re.findall(r'\d+', md_field)))  # Count digits
    mismatches = len(re.findall(r'[A-Za-z]', md_field))  # Count letters
    identity = 100 * (matches / (matches + mismatches))
    return round(identity, 2)

def sam_info_extraction(sam_file):
    sam_data = {}
    sample_name = sam_file.split("\\")[-1]
    
    with open(sam_file, 'r') as tmp:
        sam_reader = Reader(tmp)
        sam = next(sam_reader)
        lines = list(tmp)
        
        for line in tqdm(lines, desc=f"Processing {sample_name}"):
            if not line.startswith('@'):
                fields = line.strip().split('\t')
                # qname = sam.qname
                qname = fields[0]
                # ref_name = sam.rname
                ref_name = fields[2]

                if ref_name == "*":
                    continue
                
                start_position = int(fields[3])
                end_position = start_position + len(sam.seq)
                alignment_length = end_position - start_position
                md_field = None

                for field in fields[11:]:
                    if field.startswith("MD:Z:"):
                        md_field = field[5:]
                        break

                if md_field is not None:
                    percent_identity = calculate_percent_identity(md_field)
                    if ref_name not in sam_data:
                        sam_data[ref_name] = []
                    sam_data[ref_name].append({
                        'Sample': qname,
                        'Percent_Identity': percent_identity,
                        'Start_Position': start_position,
                        'End_Position': end_position,
                        'Alignment_Length': alignment_length
                    })
    return sam_data

def psl_info_extraction(psl_file):
    header_lines = 5
    psl_data = {}
    file_name = str(psl_file)
    with open(psl_file, 'r') as tmp:
        for i in range(header_lines):
            next(tmp)
        for line in tqdm(tmp, desc=f"Processing {file_name}"):
            rows = line.strip().split()
            qname = rows[9]
            ref_name = rows[13]
            matches = int(rows[0])
            mismatches = int(rows[1])
            identity = round(100 * (matches / (matches + mismatches))) # calculate percent identity
            start_pos = rows[15]
            end_pos = rows[16]
            alignment_length = int(end_pos) - int(start_pos)
            
            if ref_name not in psl_data:
                psl_data[ref_name] = []
            psl_data[ref_name].append({
                'Sample': qname,
                'Percent_Identity': identity,
                'Start_Position': start_pos,
                'End_Position': end_pos,
                'Alignment_Length': alignment_length
            })
    return psl_data
        
def frp(data_dict, output_path):
    all_data = data_dict
    output_folder = output_path
    for ref_name, data in all_data.items():
        print(f"\nPlotting FRP for {ref_name}")
        
        plt.figure(figsize=(12, 6))
            
        positions = [entry['Start_Position'] for entry in data]
        percent_identities = [entry['Percent_Identity'] for entry in data]

        plt.scatter(positions, percent_identities, s=1)
        plt.xlabel("Alignment Position")
        plt.xscale("linear")
        plt.ylabel("Percent Identity")
        plt.title(f"Fragment Recruitment Plot for {ref_name}")
        plot_file_name = os.path.join(output_folder, f"{ref_name}_{file_type}_FRP.png")
        plt.savefig(plot_file_name, dpi=300)
        plt.close()
        print(f"FRP plot saved to {plot_file_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate FRPs from sam and psl files")
    parser.add_argument("--sam_dir", type=str, help="Path to the directory containing SAM files.")
    parser.add_argument("--psl_dir", type=str, help="Path to the directory containing PSL files.")
    parser.add_argument("--output_dir", type=str, help="Path to the output directory for FRP plots.")
    args = parser.parse_args()
    
    sam_files_dir = args.sam_dir
    psl_files_dir = args.psl_dir
    output_folder = args.output_dir
    
    if sam_files_dir:
        all_sam_data = {}
        file_type = "sam"
        for file in os.listdir(sam_files_dir):
            if file.endswith(".sam"):
                file_path = os.path.join(sam_files_dir, file)
                data = sam_info_extraction(file_path)
                for ref_name, entries in data.items():
                    if ref_name not in all_sam_data:
                        all_sam_data[ref_name] = []
                    all_sam_data[ref_name].extend(entries)
        frp(all_sam_data, output_folder)
    
    elif psl_files_dir:
        all_psl_data = {}
        file_type = "psl"
        for dir in os.listdir(psl_files_dir):
            sub_dir = os.path.join(psl_files_dir, dir)
            for f in os.listdir(sub_dir):
                if f.endswith('.psl'):
                    file_path = os.path.join(sub_dir, f)
                    data = psl_info_extraction(file_path)
                    for ref_name, entries in data.items():
                        if ref_name not in all_psl_data:
                            all_psl_data[ref_name] = []
                        all_psl_data[ref_name].extend(entries)
            frp(all_psl_data, output_folder)