#!/usr/bin/env python3
    
import subprocess

# file_paths = [
#     './bowtie2_docker_image/refgenomes/tax_id_187327/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_334390/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_334413/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_391904/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_431946/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_457412/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_469604/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_557433/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_565042/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_585535/genomic.fna',
#     './bowtie2_docker_image/refgenomes/tax_id_657309/genomic.fna'
# ]

file_paths = """./bowtie2_docker_image/refgenomes/tax_id_187327/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_334390/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_334413/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_391904/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_431946/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_457412/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_469604/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_557433/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_565042/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_585535/genomic.fna,
                ./bowtie2_docker_image/refgenomes/tax_id_657309/genomic.fna"""


bowtie2_build = [
    'bowtie2-build',
    '-f', file_paths, # input file is fasta
    'ref_index'
    ]

# run using subprocess
subprocess.run(bowtie2_build)
print("Reference index built.")