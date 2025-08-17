### Parsing output from Paml
### indel_parsePaml.2.py
### v2 2024/10/16
### Parsing result from Paml output with large and small gap length
### v1 2024/10/02
### Parsing result from Paml output 
### @author:Chen Yu Chi

import os,sys
import glob
import re
import csv
import pandas as pd
import numpy as np

dir_work = '/work/myucchen/indel/indel.03.02/'

# import dataframe
fileOut = glob.glob(dir_work + '*/*/out/*/*.out')

# Regular expression pattern to extract species and their corresponding values
pat_dnds = r'([a-zA-Z_]+\d+): ([0-9.]+)'
pat_w = r'([a-zA-Z_]+\d+)\s*#\s*([0-9.]+)'

# Initialize a list to hold all data
df_m1m2 = []
df_m0 = []

# Function to extract tree values
def extract_tree_values(lines, line_index, indel, model, updwn, np, lnL, tree_label, pattern):
    if line_index + 1 < len(lines):
        tree_line = lines[line_index + 1].strip()
        species_values = re.findall(pattern, tree_line)
        for species, value in species_values:
            df_m1m2.append({
                'indel': indel,
                'model': model,
                'updwn': updwn,
                'np': np,
                'lnL': lnL,
                'species': species,
                'value': float(value),
                'type': tree_label
            })

# m1, m2 dataset
for file in fileOut:
    parsePath = re.search(r"(del\d*|ins\d*|all\d*)/(m0|m1|m2)/out/(dwn\d+|up\d+|pos\d+)/", file)
    indel, model, updwn = parsePath.groups()  # Unpack the matched groups
    if (model == 'm1') or (model == 'm2'):
        with open(file, 'r+') as outIn:
            lines = outIn.readlines()
            current_data = {}
            for i, line in enumerate(lines):
                # Extract lnL and np
                if 'lnL' in line:
                    line_lnL = re.search(r"np:\s*(\d+)\):\s*(-?\d+\.\d+)", line)
                    current_data['np'] = int(line_lnL.group(1))
                    current_data['lnL'] = float(line_lnL.group(2))
                # Check for dS, dN, and w tree lines
                for tree_label in ['dS', 'dN']:
                    if f'{tree_label} tree:' in line:
                        extract_tree_values(lines, i, indel, model, updwn, current_data.get('np'), current_data.get('lnL'), tree_label, pat_dnds)
                if 'w ratios' in line:
                    extract_tree_values(lines, i, indel, model, updwn, current_data.get('np'), current_data.get('lnL'), 'w', pat_w)

for file in fileOut:
    parsePath = re.search(r"(del\d*|ins\d*|all\d*)/(m0|m1|m2)/out/(dwn\d+|up\d+|pos\d+)/", file)
    indel, model, updwn = parsePath.groups()  # Unpack the matched groups
    if (model == 'm0'):
        with open(file, 'r+') as outIn:
            lines = outIn.readlines()
            for i, line in enumerate(lines):
                # Extract lnL and np
                if 'lnL' in line:
                    line_lnL = re.search(r"np:\s*(\d+)\):\s*(-?\d+\.\d+)", line)
                    np = int(line_lnL.group(1))
                    lnL = float(line_lnL.group(2))
                # Check for dS, dN, and w tree lines 
                if 'omega' in line:
                    line_omega = re.search(r'omega\s*\(dN/dS\)\s*=\s*([0-9.]+)', line)
                    w = float(line_omega.group(1))
                    df_m0.append({
                        'indel': indel,
                        'model': model,
                        'updwn': updwn,
                        'np': np,
                        'lnL': lnL,
                        'omega': w
                    })

# Convert df_m1m2 into a DataFrame
df_m1m2 = pd.DataFrame(df_m1m2)
df_m0 = pd.DataFrame(df_m0)

# export result
# dfOut.to_csv(dir_work + 'dnds.out', sep='\t', index=False)

df_m1m2_pivot = df_m1m2.pivot_table(index=['indel', 'model', 'updwn', 'np', 'lnL', 'type'], columns='species', values='value').reset_index()
# Optionally, flatten the MultiIndex columns
df_m1m2_pivot.columns.name = None  # Remove the name of the columns
df_m1m2_pivot.columns = df_m1m2_pivot.columns.to_flat_index()  # Flatten MultiIndex to a single index
df_m1m2_pivot.columns = [str(col) for col in df_m1m2_pivot.columns] 
df_m1m2_pivot.to_csv(dir_work + 'dnds.m1m2.out', sep='\t', index=False)
df_m0.to_csv(dir_work + 'dnds.m0.out', sep='\t', index=False)