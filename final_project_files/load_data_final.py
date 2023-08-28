#!/usr/local/bin/python3

import sqlite3 as sql
import pandas as pd

# create database in the base final project folder and connect
conn = sql.connect('C:/Users/tuann/OneDrive - Johns Hopkins/School/AS.410.712 (Advanced Practical Computer Concepts for Bioinformatics/Final Project/base/final_3.db')
curs = conn.cursor()

# set this option so I could see the columns in the console when printing
pd.set_option('display.max_columns', None)

# I put them in a list and just iterated over the list
queries = [
    """
    CREATE TABLE drug (
        drug_id varchar(200),
        drug_name varchar(200)
    );
    """,
    """
    CREATE TABLE drug_name_alt (
        drug_id varchar(200),
        drug_name_alt varchar(200)
    );
    """,
    """
    CREATE TABLE cell_line (
        cell_line_id varchar(200),
        cell_line_name varchar(200)
    );
    """,
    """
    CREATE TABLE pathway (
        path_id varchar(200),
        drug_target_pathway varchar(200)
    );
    """,
    """
    CREATE TABLE target (
        target_id varchar(200),
        drug_target varchar(200)
    );
    """,
    """
    CREATE TABLE ic50 (
        entry_id varchar(200),
        ic50 varchar(200)
    );
    """,
    """
    CREATE TABLE keys (
        entry_id varchar(200),
        drug_id varchar(200),
        cell_line_id varchar(200),
        path_id varchar(200),
        target_id varchar(200)
    );
    """
]

for query in queries:
    curs.execute(query)
    
conn.commit()

# read in the data from CSV/xlsx
# GDSC 1: cell line name, drug name, drug target, pathway, IC50
df_GDSC1 = pd.read_excel("C:\\Users\\tuann\\OneDrive - Johns Hopkins\\School\\AS.410.712 (Advanced Practical Computer Concepts for Bioinformatics\\Final Project\\base\\data\\GDSC1_fitted_dose_response_24Jul22.xlsx", usecols = ['CELL_LINE_NAME','DRUG_NAME','PUTATIVE_TARGET','PATHWAY_NAME','IC50'])
df_drug_names = pd.read_csv("C:\\Users\\tuann\\OneDrive - Johns Hopkins\\School\\AS.410.712 (Advanced Practical Computer Concepts for Bioinformatics\\Final Project\\base\\data\\screened_compounds_rel_8.4.csv", usecols = ['DRUG_NAME','SYNONYMS'])

############################
############################

# convert drug names in both dataframes to lowercase
df_drug_names['DRUG_NAME'] = df_drug_names['DRUG_NAME'].str.lower()
df_GDSC1['DRUG_NAME'] = df_GDSC1['DRUG_NAME'].str.lower()

# get unique drug names from df_GDSC1
unique_drug_names = df_GDSC1['DRUG_NAME'].unique()
print(unique_drug_names)

# filter df_drug_names to keep only the rows with drug_name in unique_drug_names
df_drug_names_filtered = df_drug_names[df_drug_names['DRUG_NAME'].isin(unique_drug_names)]
print(df_drug_names_filtered)

# make IDs for each drug name and store in the dictionary
drug_id_dict = {}

for idx, drug_name in enumerate(unique_drug_names):
    drug_id = f'{idx+1:06d}'
    drug_id_dict[drug_name] = drug_id

# add ids for the drugs that are in gdsc1
df_GDSC1['drug_id'] = df_GDSC1['DRUG_NAME'].apply(lambda x: drug_id_dict[x])    
# print(df_GDSC1)
df_drug_names_filtered['drug_id'] = df_drug_names_filtered['DRUG_NAME'].apply(lambda x: drug_id_dict[x])
print(df_drug_names_filtered)

# split the synonyms into separate rows, one for each synonym
df_drug_names_filtered['SYNONYMS'] = df_drug_names_filtered['SYNONYMS'].str.split(', ')
df_drug_names_filtered = df_drug_names_filtered.explode('SYNONYMS')
print(df_drug_names_filtered)

###### adding ids for other variables #######
df_final = df_GDSC1
df_final['cell_line_id'] = (df_final.groupby('CELL_LINE_NAME').ngroup()+1).apply(lambda x: f"{x:06}")                            
df_final['path_id'] = (df_final.groupby('PATHWAY_NAME').ngroup()+1).apply(lambda x: f"{x:03}")
df_final['target_id'] = (df_final['PUTATIVE_TARGET'].fillna('NaN').groupby(df_final['PUTATIVE_TARGET'].fillna('NaN')).ngroup() + 1).apply(lambda x: f"{x:03}")
df_final['entry_id'] = (df_final.index + 1).map(lambda x: f"{x:06}")
print(df_final)

###### queries for inserting values into db ######
# insert into drug_name_alt
for row in df_drug_names_filtered.itertuples(index=False):
    conn.execute("INSERT INTO drug_name_alt (drug_id, drug_name_alt) VALUES (?, ?)", (row.drug_id, row.SYNONYMS))

# insert into drug_name
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO drug (drug_id, drug_name) VALUES (?,?)", (row.drug_id, row.DRUG_NAME))

# insert into cell_line
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO cell_line (cell_line_id, cell_line_name) VALUES (?, ?)", (row.cell_line_id, row.CELL_LINE_NAME))
    
# insert into pathway
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO pathway (path_id, drug_target_pathway) VALUES (?,?)", (row.path_id, row.PATHWAY_NAME))
    
# insert into target
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO target (target_id, drug_target) VALUES (?, ?)", (row.target_id, row.PUTATIVE_TARGET))
    
# insert into ic50
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO ic50 (entry_id, ic50) VALUES (?, ?)", (row.entry_id, row.IC50))
    
# insert ids into keys
for row in df_final.itertuples(index=False):
    conn.execute("INSERT INTO keys (entry_id, drug_id, cell_line_id, path_id, target_id) VALUES (?,?,?,?,?)", (row.entry_id, row.drug_id, row.cell_line_id, row.path_id, row.target_id))
  
##### delete all the duplicates from cell_line, drug, pathway, target tables #####
queries_delete = [
    """
    DELETE FROM drug
    WHERE rowid NOT IN
        (SELECT min(rowid)
        FROM drug
        GROUP BY drug_name);
    """,
    """
    DELETE FROM cell_line
    WHERE rowid NOT IN
        (SELECT min(rowid)
        FROM cell_line
        GROUP BY cell_line_name);
    """,
    """
    DELETE FROM pathway
    WHERE rowid NOT IN
        (SELECT min(rowid)
        FROM pathway
        GROUP BY drug_target_pathway);
    """,
    """
    DELETE FROM target
    WHERE rowid NOT IN
        (SELECT min(rowid)
        FROM target
        GROUP BY drug_target);
    """
]

for qry in queries_delete:
    conn.execute(qry)

conn.commit()
curs.close()
conn.close()