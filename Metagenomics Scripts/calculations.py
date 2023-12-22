import os
import re
import matplotlib.pyplot as plt
import pandas as pd
from simplesam import Reader
from tqdm import tqdm
from Bio import SeqIO

sam_files_dir = "/bowtie2_local_output"
psl_files_dir = "/blat_output"

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

def average_rank(rankings):
    return sum(rankings) / len(rankings)

all_sam_data = {}

for file in os.listdir(sam_files_dir):
    if file.endswith(".sam"):
        file_path = os.path.join(sam_files_dir, file)
        data = sam_info_extraction(file_path)
        for ref_name, entries in data.items():
            if ref_name not in all_sam_data:
                all_sam_data[ref_name] = []
            all_sam_data[ref_name].extend(entries)

print(all_sam_data.keys())

sam_alignment_stats = []

for ref_name, data in all_sam_data.items():
    print(f"Calculating stats for {ref_name}")
    num_alignments = len(data) # how many entries are there for this ref
    avg_percent_identity = sum([entry['Percent_Identity'] for entry in data]) / num_alignments
    avg_alignment_length = sum([entry['Alignment_Length'] for entry in data]) / num_alignments
    
    sam_alignment_stats.append(
        {
            'Ref_Name': ref_name,
            'Num_Alignments': num_alignments,
            'Avg_Percent_Identity': avg_percent_identity,
            'Avg_Alignment_Length': avg_alignment_length
        })
    
    print(f"Number of alignments: {num_alignments}")
    print(f"Average percent identity: {avg_percent_identity}")
    print(f"Average alignment length: {avg_alignment_length}")

sam_alignment_stats_df = pd.DataFrame(sam_alignment_stats)
sam_alignment_stats_df.sort_values(by=['Num_Alignments'], ascending=False, inplace=True)
sam_alignment_stats_df.reset_index(drop=True, inplace=True)
sam_alignment_stats_df.index += 1
print(sam_alignment_stats_df)
# sorted by number of alignments
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1         NC_011593.1          370589             89.125080            351.116649
# 2         NC_010609.1          238341             87.584522            315.868457
# 3         NC_010610.1          220222             90.401439            303.807131
# 4         NC_010376.1          216994             89.258834            347.868669
# 5         NC_013654.1          157384             90.212987            371.318584
# 6         NC_015067.1          138008             89.136425            345.754703
# 7       NZ_LR698962.1           96997             88.284194            326.786509
# 8         NC_017360.1           29839             89.186534            333.010490
# 9   NZ_ACII02000003.1           29205             87.113907            352.680568
# 10  NZ_ACII02000002.1           27258             86.745230            349.069594
# 11        NC_022196.1           20975             92.537926            296.868796
# 12  NZ_ACII02000001.1            7648             87.034552            334.751438

# sorted by average percent identity
sam_alignment_stats_df.sort_values(by=['Avg_Percent_Identity'], ascending=False, inplace=True)
sam_alignment_stats_df.reset_index(drop=True, inplace=True)
sam_alignment_stats_df.index += 1
print(sam_alignment_stats_df)
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1         NC_022196.1           20975             92.537926            296.868796
# 2         NC_010610.1          220222             90.401439            303.807131
# 3         NC_013654.1          157384             90.212987            371.318584
# 4         NC_010376.1          216994             89.258834            347.868669
# 5         NC_017360.1           29839             89.186534            333.010490
# 6         NC_015067.1          138008             89.136425            345.754703
# 7         NC_011593.1          370589             89.125080            351.116649
# 8       NZ_LR698962.1           96997             88.284194            326.786509
# 9         NC_010609.1          238341             87.584522            315.868457
# 10  NZ_ACII02000003.1           29205             87.113907            352.680568
# 11  NZ_ACII02000001.1            7648             87.034552            334.751438
# 12  NZ_ACII02000002.1           27258             86.745230            349.069594

# sorted by average alignment length
sam_alignment_stats_df.sort_values(by=['Avg_Alignment_Length'], ascending=False, inplace=True)
sam_alignment_stats_df.reset_index(drop=True, inplace=True)
sam_alignment_stats_df.index += 1
print(sam_alignment_stats_df)
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1         NC_013654.1          157384             90.212987            371.318584
# 2   NZ_ACII02000003.1           29205             87.113907            352.680568
# 3         NC_011593.1          370589             89.125080            351.116649
# 4   NZ_ACII02000002.1           27258             86.745230            349.069594
# 5         NC_010376.1          216994             89.258834            347.868669
# 6         NC_015067.1          138008             89.136425            345.754703
# 7   NZ_ACII02000001.1            7648             87.034552            334.751438
# 8         NC_017360.1           29839             89.186534            333.010490
# 9       NZ_LR698962.1           96997             88.284194            326.786509
# 10        NC_010609.1          238341             87.584522            315.868457
# 11        NC_010610.1          220222             90.401439            303.807131
# 12        NC_022196.1           20975             92.537926            296.868796

# summary stats for NC_010610.1
data_for_target_ref = all_sam_data['NC_010610.1']
target_ref_df = pd.DataFrame(data_for_target_ref)

# target_ref_df['Percent_Identity'].describe()
target_ref_df['Alignment_Length'].describe()
# NC_010610.1 alignment length summary stats
# count    220222.000000
# mean        303.807131
# std         202.781976
# min          17.000000
# 25%          27.000000
# 50%         368.000000
# 75%         489.000000
# max         544.000000

sam_alignment_stats_df['Rank_Num_Alignments'] = sam_alignment_stats_df['Num_Alignments'].rank(ascending=False)
sam_alignment_stats_df['Rank_Avg_Percent_Identity'] = sam_alignment_stats_df['Avg_Percent_Identity'].rank(ascending=False)
sam_alignment_stats_df['Rank_Avg_Alignment_Length'] = sam_alignment_stats_df['Avg_Alignment_Length'].rank(ascending=False)
sam_alignment_stats_df['Avg_Rank'] = sam_alignment_stats_df[['Rank_Num_Alignments', 'Rank_Avg_Percent_Identity', 'Rank_Avg_Alignment_Length']].apply(average_rank, axis=1)
print(sam_alignment_stats_df[['Ref_Name', 'Avg_Rank']].sort_values(by='Avg_Rank'))
# sorted by average rank
#              Ref_Name   Avg_Rank
# 5         NC_013654.1   3.000000
# 0         NC_011593.1   3.666667
# 2         NC_010376.1   4.333333
# 8         NC_010610.1   5.333333
# 1         NC_015067.1   6.000000
# 3         NC_010609.1   7.000000
# 4   NZ_ACII02000003.1   7.000000
# 10        NC_017360.1   7.000000
# 7       NZ_LR698962.1   8.000000
# 9         NC_022196.1   8.000000
# 6   NZ_ACII02000002.1   8.666667
# 11  NZ_ACII02000001.1  10.000000

sam_alignment_stats_df.sort_values(by='Ref_Name')
sam_alignment_stats_df

all_psl_data = {}

for folder in os.listdir(psl_files_dir):
    # print(f"Processing {folder}")
    sub_dir = os.path.join(psl_files_dir, folder)
    # print(sub_dir)
    for f in os.listdir(sub_dir):
        if f.endswith('.psl'):
            file_path = os.path.join(sub_dir, f)
            data = psl_info_extraction(file_path)
            for ref_name, entries in data.items():
                if ref_name not in all_psl_data:
                    all_psl_data[ref_name] = []
                all_psl_data[ref_name].extend(entries)

print(all_psl_data.keys())

psl_alignment_stats = []

for ref_name, data in all_psl_data.items():
    print(f"Calculating stats for {ref_name}")
    num_alignments = len(data) # how many entries are there for this ref
    avg_percent_identity = sum([entry['Percent_Identity'] for entry in data]) / num_alignments
    avg_alignment_length = sum([entry['Alignment_Length'] for entry in data]) / num_alignments
    
    psl_alignment_stats.append(
        {
            'Ref_Name': ref_name,
            'Num_Alignments': num_alignments,
            'Avg_Percent_Identity': avg_percent_identity,
            'Avg_Alignment_Length': avg_alignment_length
        })
    
    print(f"Number of alignments: {num_alignments}")
    print(f"Average percent identity: {avg_percent_identity}")
    print(f"Average alignment length: {avg_alignment_length}")

psl_alignment_stats_df = pd.DataFrame(psl_alignment_stats)
psl_alignment_stats_df.sort_values(by=['Num_Alignments'], ascending=False, inplace=True)
psl_alignment_stats_df.reset_index(drop=True, inplace=True)
psl_alignment_stats_df.index += 1
print(psl_alignment_stats_df)
# sorted by number of alignments
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1         NC_013654.1         9813292             88.458128            474.599605
# 2         NC_010609.1         8246223             88.976295            522.084593
# 3         NC_022196.1         7102334             88.366525            448.443686
# 4         NC_015067.1         5547239             88.926861            535.967907
# 5         NC_011593.1         5539710             89.017243            538.179225
# 6   NZ_ACII02000003.1         4186540             88.264326           1080.864151
# 7   NZ_ACII02000002.1         4155084             88.464868            617.032779
# 8         NC_017360.1         2428391             89.976393            357.670278
# 9   NZ_ACII02000001.1         1376173             88.129484            373.645851
# 10        NC_010610.1         1262748             93.045654            371.216551
# 11        NC_010376.1         1195342             94.061447            358.535053
# 12      NZ_LR698962.1          458513             93.756469            311.463029

psl_alignment_stats_df.sort_values(by=['Avg_Percent_Identity'], ascending=False, inplace=True)
psl_alignment_stats_df.reset_index(drop=True, inplace=True)
psl_alignment_stats_df.index += 1
print(psl_alignment_stats_df)
# sorted by average percent identity
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1         NC_010376.1         1195342             94.061447            358.535053
# 2       NZ_LR698962.1          458513             93.756469            311.463029
# 3         NC_010610.1         1262748             93.045654            371.216551
# 4         NC_017360.1         2428391             89.976393            357.670278
# 5         NC_011593.1         5539710             89.017243            538.179225
# 6         NC_010609.1         8246223             88.976295            522.084593
# 7         NC_015067.1         5547239             88.926861            535.967907
# 8   NZ_ACII02000002.1         4155084             88.464868            617.032779
# 9         NC_013654.1         9813292             88.458128            474.599605
# 10        NC_022196.1         7102334             88.366525            448.443686
# 11  NZ_ACII02000003.1         4186540             88.264326           1080.864151
# 12  NZ_ACII02000001.1         1376173             88.129484            373.645851

psl_alignment_stats_df.sort_values(by=['Avg_Alignment_Length'], ascending=False, inplace=True)
psl_alignment_stats_df.reset_index(drop=True, inplace=True)
psl_alignment_stats_df.index += 1
print(psl_alignment_stats_df)
# sorted by average alignment length
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 1   NZ_ACII02000003.1         4186540             88.264326           1080.864151
# 2   NZ_ACII02000002.1         4155084             88.464868            617.032779
# 3         NC_011593.1         5539710             89.017243            538.179225
# 4         NC_015067.1         5547239             88.926861            535.967907
# 5         NC_010609.1         8246223             88.976295            522.084593
# 6         NC_013654.1         9813292             88.458128            474.599605
# 7         NC_022196.1         7102334             88.366525            448.443686
# 8   NZ_ACII02000001.1         1376173             88.129484            373.645851
# 9         NC_010610.1         1262748             93.045654            371.216551
# 10        NC_010376.1         1195342             94.061447            358.535053
# 11        NC_017360.1         2428391             89.976393            357.670278
# 12      NZ_LR698962.1          458513             93.756469            311.463029

psl_alignment_stats_df['Rank_Num_Alignments'] = psl_alignment_stats_df['Num_Alignments'].rank(ascending=False)
psl_alignment_stats_df['Rank_Avg_Percent_Identity'] = psl_alignment_stats_df['Avg_Percent_Identity'].rank(ascending=False)
psl_alignment_stats_df['Rank_Avg_Alignment_Length'] = psl_alignment_stats_df['Avg_Alignment_Length'].rank(ascending=False)
psl_alignment_stats_df['Avg_Rank'] = psl_alignment_stats_df[['Rank_Num_Alignments', 'Rank_Avg_Percent_Identity', 'Rank_Avg_Alignment_Length']].apply(average_rank, axis=1)
print(psl_alignment_stats_df[['Ref_Name', 'Avg_Rank']].sort_values(by='Avg_Rank'))
# sorted by average rank
#              Ref_Name  Avg_Rank
# 2         NC_010609.1  4.333333
# 5         NC_011593.1  4.333333
# 4         NC_015067.1  5.000000
# 1         NC_013654.1  5.333333
# 7   NZ_ACII02000002.1  5.666667
# 6   NZ_ACII02000003.1  6.000000
# 3         NC_022196.1  6.666667
# 10        NC_010610.1  7.333333
# 11        NC_010376.1  7.333333
# 8         NC_017360.1  7.666667
# 12      NZ_LR698962.1  8.666667
# 9   NZ_ACII02000001.1  9.666667

psl_alignment_stats_df.sort_values(by='Ref_Name')

# sam
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length
# 5         NC_010376.1          216994             89.258834            347.868669
# 10        NC_010609.1          238341             87.584522            315.868457
# 11        NC_010610.1          220222             90.401439            303.807131
# 3         NC_011593.1          370589             89.125080            351.116649
# 1         NC_013654.1          157384             90.212987            371.318584
# 6         NC_015067.1          138008             89.136425            345.754703
# 8         NC_017360.1           29839             89.186534            333.010490
# 12        NC_022196.1           20975             92.537926            296.868796
# 7   NZ_ACII02000001.1            7648             87.034552            334.751438
# 4   NZ_ACII02000002.1           27258             86.745230            349.069594
# 2   NZ_ACII02000003.1           29205             87.113907            352.680568
# 9       NZ_LR698962.1           96997             88.284194            326.786509

# psl
#              Ref_Name  Num_Alignments  Avg_Percent_Identity  Avg_Alignment_Length  Rank_Num_Alignments  Rank_Avg_Percent_Identity  Rank_Avg_Alignment_Length  Avg_Rank
# 11        NC_010376.1         1195342             94.061447            358.535053                 11.0                        1.0                       10.0  7.333333    
# 2         NC_010609.1         8246223             88.976295            522.084593                  2.0                        6.0                        5.0  4.333333    
# 10        NC_010610.1         1262748             93.045654            371.216551                 10.0                        3.0                        9.0  7.333333    
# 5         NC_011593.1         5539710             89.017243            538.179225                  5.0                        5.0                        3.0  4.333333    
# 1         NC_013654.1         9813292             88.458128            474.599605                  1.0                        9.0                        6.0  5.333333    
# 4         NC_015067.1         5547239             88.926861            535.967907                  4.0                        7.0                        4.0  5.000000    
# 8         NC_017360.1         2428391             89.976393            357.670278                  8.0                        4.0                       11.0  7.666667    
# 3         NC_022196.1         7102334             88.366525            448.443686                  3.0                       10.0                        7.0  6.666667    
# 9   NZ_ACII02000001.1         1376173             88.129484            373.645851                  9.0                       12.0                        8.0  9.666667    
# 7   NZ_ACII02000002.1         4155084             88.464868            617.032779                  7.0                        8.0                        2.0  5.666667    
# 6   NZ_ACII02000003.1         4186540             88.264326           1080.864151                  6.0                       11.0                        1.0  6.000000    
# 12      NZ_LR698962.1          458513             93.756469            311.463029                 12.0                        2.0                       12.0  8.666667 