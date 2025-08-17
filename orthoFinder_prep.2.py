### Organizer orthofinder result
### orthoFinder_prep.2.py
### v2 2024/05/28
### Organize the result from orthofinder with agat processed (only leave the longest isoform)
### v1 2024/05/10
### Organize the result from orthofinder to make list for extract sequences and relabel sequences
### @author:Chen Yu Chi

import os,sys
import glob
import re
import csv
import pandas as pd
# pyjanitor
import janitor
from Bio import SeqIO
from Bio import AlignIO

dir_in = '/work/myucchen/indel/indel.02.02/'
dir_ortho = dir_in + 'orthofinder/Results_Aug16/'
file_ortho = dir_ortho + 'Orthogroups/Orthogroups.tsv'
dir_cds = dir_in + 'cds/'
dir_aa = dir_in + 'aa/'
dir_out = '/work/myucchen/indel/indel.02.03/'
lst_sp = '/work/myucchen/indel/indel.02.01/species.lst'
fa_ext = 'fasta'

os.system('mkdir -p ' + dir_out + 'cds/')
os.system('mkdir -p ' + dir_out + 'aa/')

# import dataframe
df_ortho = pd.read_csv(file_ortho, sep = '\t', header = 0)
# remove rows with column of NA and remove space
df_ortho = df_ortho.dropna() # any species missing in the orthogroups
df_ortho = df_ortho.replace(' ', '', regex =True)
# melt the dataframe
df_ortho = df_ortho.pivot_longer(index = 'Orthogroup', names_to = ('species'), values_to = ('rna'))
df_ortho['rna'] = df_ortho['rna'].str.split(',')
df_ortho = df_ortho.explode('rna').reset_index(drop=True)
df_ortho['rna'] = df_ortho['rna'].replace('rna-', '', regex =True)

df_ortho['abbr'] = df_ortho['species'].str.split('_', expand=True)[0].str[0:1] + df_ortho['species'].str.split('_', expand=True)[1].str[0:2]
df_ortho['abbr'] = df_ortho['abbr'].str.lower()
#df_ortho.to_csv(dir_out + 'ortho.lst', sep = '\t', index = None)

# Parse agat seq info
lst_fa = glob.glob(dir_cds + '*.' + fa_ext)

# patterns
pat_rna = r'rna-([^ ]+) '
pat_gen = r'gene=gene-([^ ]+) '
pat_seqID = r'seq_id=([^ ]+) '
pat_type = r'type=([^ ]+)'

# compile all patterns
patterns = {
    'rna': pat_rna,
    'gene': pat_gen,
    'seqID': pat_seqID,
    'type': pat_type,
}

# Functions for extracting info
def extract_info(description, patterns):
    info = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, description)
        info[key] = match.group(1) if match else None
    return info

df_agat = []
for file in lst_fa:
    species = re.sub('^.*/', '', file)
    species = re.sub('.fasta', '', species)
    print(species)
    with open(file, 'r') as inSeq:
        for seq in SeqIO.parse(inSeq, 'fasta'):
            info = extract_info(seq.description, patterns)
            #print(info)
            info['species'] = species
            df_agat.append(info)

df_agat = pd.DataFrame(df_agat)
df_ortho = pd.merge(df_ortho, df_agat, left_on = ["rna", "species"], right_on = ["rna", "species"], how = "left")
df_ortho['file_name'] = df_ortho['Orthogroup'].str.strip('OG').astype(int).astype(str)
df_ortho['sn'] = df_ortho.groupby(['Orthogroup', 'species']).cumcount() + 1
df_ortho['seq_name'] = df_ortho['abbr'] + '_' + df_ortho['sn'].astype(str)
df_ortho.to_csv(dir_out + 'seqName.lst', sep = '\t', index = None)

# Prepare orthogroups cds files
lst_fa = glob.glob(dir_ortho + 'Orthogroup_Sequences/*.fa')

for fa in lst_fa:
    id_file = re.sub('^.*OG', 'OG', fa)
    id_file = re.sub('.fa', '', id_file)
    if id_file not in df_ortho['Orthogroup'].unique(): continue
    file_out = df_ortho[df_ortho['Orthogroup'] == id_file]['file_name'].values[0]
    #print(file_out)
    #print(id_file)
    with open(fa, 'r') as faIn, open(dir_out + 'cds/' + file_out + '.' + fa_ext, 'a+') as faOut:
        for seq in SeqIO.parse(faIn, 'fasta'):
            seq_id = re.sub('rna-', '', seq.id)
            seq.id = df_ortho[df_ortho['rna'] == seq_id]['seq_name'].values[0]
            #print(seq.id))
            faOut.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')

# Prepare orthogroups aa files
lst_fa = glob.glob(dir_aa + '*.fa')
df_ortho['rna_id'] = df_ortho['rna'].str.replace(r'.{4}$', '', regex=True)
df_sp = pd.read_csv(lst_sp, sep = '\t', header = None, names = ['sp_name', 'id', 'abbr'])


for fa in lst_fa:
    fa_name = re.sub('^.*/', '', fa)
    fa_name = fa_name[:15]
    sp =  df_sp[df_sp['id'] == fa_name]['abbr'].values[0]
    print(fa_name + ' ' + sp)
    with open(fa) as faIn:
        for seq in SeqIO.parse(faIn, 'fasta'):
            seq_id = re.sub('rna-', '', seq.id)
            seq_id2 = seq_id + '_' + sp
            #print(seq_id2)
            if seq_id2 not in df_ortho['rna'].unique(): continue
            file_out = df_ortho[df_ortho['rna'] == seq_id2]['file_name'].values[0]
            seq.id = df_ortho[df_ortho['rna'] == seq_id2]['seq_name'].values[0]
            #print(seq.id)
            with open(dir_out + 'aa/' + file_out + '.' + fa_ext, 'a+') as faOut:
                faOut.write('>' + str(seq.id) + '\n' + str(seq.seq) + '\n')
