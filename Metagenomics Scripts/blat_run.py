#!/usr/bin/env python3

import os
import subprocess
import re
from tqdm import tqdm

file_paths = ["/home/tuann99/metagenomics/midterm/refgenomes/tax_id_187327/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_334390/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_334413/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_391904/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_431946/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_457412/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_469604/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_557433/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_565042/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_585535/genomic.fna",
              "/home/tuann99/metagenomics/midterm/refgenomes/tax_id_657309/genomic.fna"]
data = "/home/tuann99/metagenomics/midterm/data/"
output = "/home/tuann99/metagenomics/midterm/blat_output/"

for path in file_paths: # aligning against each ref separately
    tax_id = re.search(r'tax_id_([0-9]+)', path).group(1) # get tax id from path
    tax_id_output_dir = os.path.join(output, tax_id) # create output dir for each tax id
    os.makedirs(tax_id_output_dir, exist_ok=True) # create dir if not exist
    tqdm.write(f'Using tax id {tax_id} as reference.')
    
    for file in tqdm(os.listdir(data)): 
        if file.endswith(".fsa"): # aligning each sample against each ref
            sample_name = re.search(r'([A-Za-z0-9]+).fsa', file).group(1) # get sample name from path
            output_path_final = os.path.join(tax_id_output_dir, f'{sample_name}_{tax_id}-as-ref_blat_output.psl') # create output path
            data_path = os.path.join(data, file) # create data path
            
            # command line arguments for blat
            blat = [
                'blat',
                '-q=dna',
                '-minIdentity=80',
                path,
                data_path,
                output_path_final
            ] 
            # use -nohead option next time to remove psl header
            
            subprocess.run(blat)
            # sys.stdout.flush()
            tqdm.write(f"Sample {sample_name} aligned.")      
