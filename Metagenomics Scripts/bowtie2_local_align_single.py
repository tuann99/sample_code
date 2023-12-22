import subprocess
import os
import time

# data_directory = './bowtie2_docker_image/data'
# output_directory = './bowtie2_docker_image/bowtie2_local_output'
data_directory = './metagenomics/midterm/data'
output_directory = './metagenomics/midterm/bowtie2_local_output'

file_counter = 0 
os.makedirs(output_directory, exist_ok=True)
start_time = time.time()

for file in os.listdir(data_directory):
    if file.endswith(".fsa"):
        file_path = os.path.join(data_directory, file)
        file_name = file.split(".")[0]
        file_counter += 1
        print(f"Processing file {file_counter}/{len(os.listdir(data_directory))}: {file_name}...")
        output_file_path = os.path.join(output_directory, f"{file_name}.sam")
        with open(output_file_path, 'w'): # for creating the file bc was getting errors without it
            pass
        
        bowtie2_local = [
            'bowtie2',
            '-f', # input file is fasta
            '-t',
            '-p', '10', # number of threads
            '-x', 'ref_index',  
            '-U', file_path,
            '--local',
            '-S', output_file_path
            ]
        
        # run using subprocess
        subprocess.run(bowtie2_local)
        print(f"Local alignment of {file_name} complete.\n")
        
end_time = time.time()
total_runtime = end_time - start_time