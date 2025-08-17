### alignment and ancestral reconstruction
### indel_aln_recon.3.py
### v3 2024/12/18
### remove overlap surrounding regions and create control
### v2 2024/08/06
### organized
### v1 2024/07/06
### align and reconstruct ancestral sequences
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

dir_work = '/work/myucchen/indel/indel.03.01/'
dir_in = '/work/myucchen/indel/indel.02.04/seq/'
dir_seq = dir_work + 'seq/'
fa_ext = 'fasta'

# make directory
os.system('mkdir -p ' + dir_work + 'bed/')
os.system('mkdir -p ' + dir_work + 'bed_cln/')
os.system('mkdir -p ' + dir_work + 'bed_ctrl/')

# copy data
#os.system('cp -R ' + dir_in + ' ' + dir_seq)

# alignment
lst_fa = glob.glob(dir_seq + '*/*_nostp.fasta')
len(lst_fa)
# 12394

# identify gaps (any gaps)
ignore_anc = ['ROOT', '1']
aln_nogap = [] # for control
aln_gap =  []
for og in lst_fa:
    dir_og = re.sub(r'[^/]*.fasta$', '', og)
    #print(dir_og)
    os.chdir(dir_og)
    id_no = re.sub('^.*/', '', og)
    id_no = re.sub(r'_.*$', '', id_no)
    print('ID: ' + id_no)
    #os.system('cat ' + id_no + '.codon ' + id_no + '_internal_ancestral_reconstruction.fas > ' +  id_no + '.' + fa_ext)
    #with open(id_no + '_nostp.' + fa_ext, 'r+') as alnIn, open(id_no + '.' + fa_ext + '_nostp.gap.bed', 'a+') as bedOut:
    #    for seq in SeqIO.parse(alnIn, 'fasta'):
    #        # skip non-interest ancestor
    #        if seq.name in ignore_anc: continue
    #        lst_gap = list(re.finditer('[-]+', str(seq.seq)))
    #        if not lst_gap: 
    #            aln_nogap.append(og + + ' ') 
    #        for region_number, match in enumerate(lst_gap, 1):
    #            bedOut.write(seq.name + '\t' + str(match.start()) + '\t' +  str(match.end()) + '\n')
    no_gaps_in_all = True  # Assume initially that all sequences have no gaps
    with open(id_no + '_nostp.' + fa_ext, 'r+') as alnIn, open(id_no + '.' + fa_ext + '_nostp.gap.bed', 'a+') as bedOut:
        for seq in SeqIO.parse(alnIn, 'fasta'):
            # Skip non-interest ancestor
            if seq.name in ignore_anc:
                continue
            if '-' in str(seq.seq): 
                no_gaps_in_all = False  
                lst_gap = list(re.finditer('[-]+', str(seq.seq)))
                for region_number, match in enumerate(lst_gap, 1):
                    bedOut.write(seq.name + '\t' + str(match.start()) + '\t' + str(match.end()) + '\n')
    if no_gaps_in_all: 
        aln_nogap.append(og) 
    else: 
        aln_gap.append(og) 

# check alinment with gaps and without gaps
len(aln_nogap)
#2134
len(aln_gap)
#10260

# output alnment list without gap for control analysis
with open(dir_work + 'aln_nogap.lst', 'w') as alnOut:
    for line in aln_nogap:
        alnOut.write(f"{line}\n")

# get the whole aln length
og_len = pd.DataFrame(columns=['OGid', 'OGlen'])
for dir_og in lst_fa:
    id_no = re.sub('^.*/', '', dir_og)
    id_no = re.sub(r'_nostp\..*$', '', id_no)
    print(id_no)
    try:
        with open(dir_og, 'r+') as alnIn:
            for seq in SeqIO.parse(alnIn, 'fasta'):
                if seq.name == 'pma_1':
                    seqLen = len(seq.seq)
                    # print(seqLen)
                    df_tmp = pd.DataFrame({'OGid':[id_no], 'OGlen': [str(seqLen)]})
                    og_len = pd.concat([og_len, df_tmp], axis = 0).reset_index(drop = True)
    except OSError:
        # You could log/print a warning here if you need.
        continue

og_len.to_csv(dir_work + 'og_len.lst', sep = '\t', index = None)

# identify insertions or deletions
# prepare bed files
os.system('cp ' + dir_seq + '*/*_nostp.gap.bed ' + dir_work + 'bed/')

# import gap
lst_bed = glob.glob(dir_work + 'bed/*.bed')
len(lst_bed)
#12394

# filter 1: indel rules
# gap in anc but not in pma: ins
# gap in pma but not in anc: del
# No gaps in the up and down stream n bp in other species 
focal_sp = ['pma_1', '2'] # focal species
buffer = 300
df_indel = pd.DataFrame(columns=['species', 'start', 'end', 'OGid'])

# Function to check for overlaps
def overlaps(start1, end1, start2, end2):
    return (start1 <= end2) and (end1 >= start2)

for og in lst_bed:
    df_bed = pd.read_csv(og, sep = '\t', names = ['species', 'start', 'end']).reset_index(drop=True)
    df_bed['start'] = pd.to_numeric(df_bed['start'])
    df_bed['end'] = pd.to_numeric(df_bed['end'])
    #print(df_bed)
    df_pma = df_bed[df_bed['species'].isin(focal_sp)].reset_index(drop=True)
    focal_gaps = []
    # gaps overlapped in focal species and overlap surrounding regions
    gaps_to_skip = set()
    for i, gap1 in df_pma.iterrows():
        for j, gap2 in df_pma.iterrows():
            if i != j :
                if overlaps(gap1['start'], gap1['end'], gap2['start'], gap2['end']):
                    gaps_to_skip.add((gap1['start'], gap1['end']))
                    gaps_to_skip.add((gap2['start'], gap2['end'])) 
    for i in range(len(df_pma)):
        gap = df_pma.iloc[i]
        start = gap['start']
        end = gap['end']
        species = gap['species']
        if (start, end) in gaps_to_skip: continue
        # Define buffer range
        buffer_start = start - buffer
        buffer_end = end + buffer
        #print(buffer_end)
        # Skip gaps if the buffer range overlaps with any gaps in any species (excluding exact matches)
        overlap_found = df_bed.apply(lambda row: (overlaps(buffer_start, start, row['start'], row['end']) or overlaps(end, buffer_end, row['start'], row['end'])) and not (start == row['start'] and end == row['end']), axis=1).any()
        if overlap_found: continue
        # If all conditions are met, add to focal gaps
        focal_gaps.append(gap)
    df_tmp = pd.DataFrame(focal_gaps)
    if df_tmp.empty: continue
    #print(df_tmp)
    id_no = re.sub('^.*/', '', og)
    id_no = re.sub('.fasta.*$', '', id_no)
    df_tmp['OGid'] = id_no
    df_indel = pd.concat([df_indel, df_tmp], axis = 0).reset_index(drop = True)
    #df_temp.to_csv(dir_work + 'bed_cln/' + id_no + 'indel.bed', sep = '\t', index = None)


"""
# Remove empty output
lst_indel = glob.glob(dir_work + 'bed_cln/*.bed')
lst_empt = []
for bed in lst_indel:
    file_size = os.path.getsize(bed)
    #print(bed, file_size)
    if file_size <= 1:
        id_no = re.sub('^.*/', '', bed)
        lst_empt.append(id_no)
        os.system('rm ' + bed)  

with open(dir_work + 'filterOut_indel.lst', 'w+') as outFile:
    outFile.write("\n".join(lst_empt))
"""

# filter 2: remove overlaping indel surrounding regions and the edge ones
# import seq length
og_len = pd.read_csv(dir_work + 'og_len.lst', sep = '\t', header = 0) 
og_len['OGid'] = og_len['OGid'].astype(str)

"""
# import bed files
lst_indel = glob.glob(dir_work + 'bed_cln/*.bed')
df_indel = pd.DataFrame()
for bed in lst_indel:
    id_no = re.sub('^.*/', '', bed)
    id_no = re.sub(r'\..*$', '', id_no)
    print(id_no)
    df_bed = pd.read_csv(bed, sep = '\t', header = 0)
    df_bed['OGid'] = id_no
    df_indel = pd.concat([df_indel, df_bed], axis = 0).reset_index(drop = True)
"""

# check if the indel is on end position
def assign_label(i, j, k):
    if i <= (0 + 3) or j >= (k - 3):
        return 'end'
    else:
        return 'mid'

# merge indel and gene length
df_indel = df_indel.merge(og_len, on = ['OGid'], how = 'left')
df_indel['OGlen'] = df_indel['OGlen'].astype(int)
df_indel['gapLen'] = df_indel['end'] - df_indel['start']
df_indel['pos'] = df_indel.apply(lambda x: assign_label(x.start, x.end, x.OGlen), axis=1)
#df_indel.to_csv(dir_work + 'indel_all.bed', sep = '\t', index = None)

# Get the indel in the middle of the genes
df_indel_mid = df_indel[df_indel['pos'] == 'mid']

# surrounding regions range
df_indel_mid["start_buff"] = df_indel_mid["start"] - buffer
df_indel_mid["end_buff"] = df_indel_mid["end"] + buffer

# remove surrounding region overlap
overlap_indices = set()
adj_reg_overlap = []

# Check overlaps within the same OGid
for OGid, group in df_indel_mid.groupby("OGid"):
    group = group.sort_values("start")  # Sort by start for efficiency
    #if OGid != '4220': continue
    for i in range(len(group)):
        for j in range(len(group)):
            row1 = group.iloc[i]
            row2 = group.iloc[j]
            if overlaps(row1["start_buff"], row1["start"], row2["end"], row2["end_buff"]) or overlaps(row2["start_buff"], row2["start"], row1["end"], row1["end_buff"]):
                overlap_indices.add(row1.name)
                overlap_indices.add(row2.name)
                adj_reg_overlap.append({
                    "OGid": OGid,
                    #"row1_index": row1.name,
                    #"row2_index": row2.name,
                    "start": row1["start"],
                    "end": row1["end"]
                })
                adj_reg_overlap.append({
                    "OGid": OGid,
                    #"row1_index": row1.name,
                    #"row2_index": row2.name,
                    "start": row2["start"],
                    "end": row2["end"]
                })


# Create DataFrames for overlapping and non-overlapping rows
adj_reg_overlap = pd.DataFrame(adj_reg_overlap)
adj_reg_overlap = adj_reg_overlap.drop_duplicates(keep = 'first')
final_indel = df_indel_mid.loc[~df_indel_mid.index.isin(overlap_indices)]

final_indel.to_csv(dir_work + 'indel_nonoverlap.bed', sep = '\t', index = None)
adj_reg_overlap.to_csv(dir_work + 'check.bed', sep = '\t', index = None)


final_indel['gapLen'] = pd.to_numeric(final_indel['gapLen'], errors='coerce')  # Converts non-numeric values to NaN
final_indel['gapLen'].describe()
final_indel['OGid'].nunique() # orthogroups with proper gaps
final_indel[final_indel['species'] == '2']['species'].count() # insertion
final_indel[final_indel['species'] == 'pma_1']['species'].count() # deletion
final_indel[final_indel['gapLen'] >= 50]['species'].count() # large indel
final_indel[final_indel['gapLen'] < 50]['species'].count() # small indel
final_indel.to_csv(dir_work + 'OGindelStat.lst', sep = '\t', index = None)


for ogid, group in final_indel.groupby('OGid'):
    # Select only the required columns
    indel = group[['species', 'start', 'end']]
    # Write to a file named after the OGid
    filename = os.path.join(dir_work, 'bed_cln/', f"{ogid}.bed")
    indel.to_csv(filename, sep='\t', index=False, header=False)


# control: orthogroup without gaps
aln_nogap = [re.sub(r'^.*/|_.*fasta', '', ogid) for ogid in aln_nogap]
indel_ctrl = pd.DataFrame({'OGid':aln_nogap})
indel_ctrl = indel_ctrl.merge(og_len, on='OGid', how='left')
indel_ctrl["OGlen"] = pd.to_numeric(indel_ctrl["OGlen"], errors="coerce")

bed_ctrl = []
for _, row in indel_ctrl.iterrows():
    oglen = row["OGlen"]
    if oglen >= buffer:
        if (oglen/3) % 2 == 0:  # Even
            start = int((oglen - buffer)*0.5)
            end = int((oglen + buffer)*0.5)
        else:  # Odd
            start = int((oglen - buffer + 3)*0.5)
            end = int((oglen + buffer + 3)*0.5)
        bed_ctrl.append({"OGid": row["OGid"], "OGlen": oglen, "start": start, "end": end})

bed_ctrl = pd.DataFrame(bed_ctrl)

# Convert result to DataFrame
lst_spe = ['pma_1', 'cco_1', 'fal_1', 'tgu_1']

for ogid, group in bed_ctrl.groupby('OGid'):
    expanded_rows = []
    for species in lst_spe:
        for _, row in group.iterrows():
            expanded_rows.append({'species': species, 'start': row['start'], 'end': row['end']})
    ctrl = pd.DataFrame(expanded_rows)
    filename = os.path.join(dir_work, 'bed_ctrl/', f"{ogid}.bed")
    ctrl.to_csv(filename, sep='\t', index=False, header=False)

