#!usr/bin/env python3

import argparse
import gzip
from Bio import SeqIO
import regex as re
import tarfile

def parse_fasta_fastq(file_path):
    try:
        is_gzipped = file_path.lower().endswith(".gz")
        is_tar_bz2 = file_path.lower().endswith(".tar.bz2")
        print("Checking if file is compressed...")

        if is_gzipped:
            print("File is gzipped.")
            handle = gzip.open(file_path, 'rt')
            f = re.sub(r'\.gz$', '', file_path, flags=re.IGNORECASE)
            print("File has been un-gzipped.")
        elif is_tar_bz2:
            print("File is tar.bz2.")
            with tarfile.open(file_path, 'r:bz2') as handle:
                f = re.sub(r'\.tar.bz2$', '', file_path, flags=re.IGNORECASE)
                print("File has been un-tarred.")
        else:
            handle = open(file_path, 'rt')
            f = file_path
            print("File is not a zipped file. Proceeding with analysis..")

        is_fasta = bool(re.search(r"\.fasta$|\.fa$", f.lower()))
        is_fastq = bool(re.search(r"\.fastq$|\.fq$", f.lower()))

        sequence_count = 0
        residue_count = 0

        if is_fasta:
            file_format = "fasta"
            print("Parsing FASTA file...")
            for record in SeqIO.parse(handle, file_format):
                sequence_count += 1
                residue_count += len(record.seq)
        elif is_fastq:
            file_format = "fastq" # only works for fastq-sanger format
            print("Parsing FASTQ file...")
            for record in SeqIO.parse(handle, file_format):
                sequence_count += 1
                residue_count += len(record.seq)
        else:
            raise ValueError("Unsupported file format. Please ensure file format is FASTA or FASTQ.")  

        print(f"File: {file_path}")
        print(f"Sequence Count: {sequence_count}")
        print(f"Residue Count: {residue_count}")
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except ValueError as e:
        print(f"Error: {e}")

def main():
    parser = argparse.ArgumentParser(description="Parse FASTA/FASTQ files and calculate sequence and residue counts.")
    parser.add_argument("file_path", help="Path to the FASTA or FASTQ file")
    args = parser.parse_args()
    parse_fasta_fastq(args.file_path)

if __name__ == "__main__":
    main()

# wget http://downloads.hmpdacc.org/data/Illumina/PHASEII/anterior_nares/SRS077085.tar.bz2 -OutFile "SRS077085"
# tar -xjf SRS077085.tar.bz2

# $ cd C:\Users\nguye620\Desktop\test_directory
# cd C:\Users\tuann\Desktop\td 
# $ gzip CAM_PROJ_SargassoSea.read_pep.fa
# $ fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.1.fastq
# $ fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.2.fastq
# $ fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.singleton.fastq

# Output:
# C:\Users\tuann\Desktop\td>python fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.1.fastq
# No zipped file found. Proceeding with analysis..
# File is FASTQ.
# File: SRS077085.denovo_duplicates_marked.trimmed.1.fastq
# Sequence Count: 58172425
# Residue Count: 5817242500

# C:\Users\tuann\Desktop\td>python fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.2.fastq 
# No zipped file found. Proceeding with analysis..
# File is FASTQ.
# File: SRS077085.denovo_duplicates_marked.trimmed.2.fastq
# Sequence Count: 58172425
# Residue Count: 5817242500

# C:\Users\tuann\Desktop\td>python fasta_fastq_summary_stats.py SRS077085.denovo_duplicates_marked.trimmed.singleton.fastq
# No zipped file found. Proceeding with analysis..
# File is FASTQ.
# File: SRS077085.denovo_duplicates_marked.trimmed.singleton.fastq
# Sequence Count: 0
# Residue Count: 0
