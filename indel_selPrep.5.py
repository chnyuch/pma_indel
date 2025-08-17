### Prepare indel files for paml
### indel_selPrep.5.py
### v5 2024/09/30
### go back to paml pipeline
### v4 2024/09/26
### concatenate the selected sequences on the same genes for hyphy
### v3 2024/08/27
### prepare sequences file for hyphy, concatenate whole sequences to one and do the hyphy selection test
### v2 2024/08/08
### improve version
### v1 2024/08/06
### prepare sequences file for paml
### @author:Chen Yu Chi

import os,sys
import glob
import re
import csv
import pandas as pd
import numpy as np
# pyjanitor
import janitor
from Bio import SeqIO, AlignIO
from Bio.Nexus import Nexus
from collections import defaultdict

dir_work = '/work/myucchen/indel/indel.03.01/'
dir_seq = dir_work + 'seq/'
dir_out = '/work/myucchen/indel/indel.03.02/inputs/'


# create nt bed for surrounding regions
# make directory
os.system('mkdir -p ' + dir_work + 'bed_surr/')


# import list
df_indel = pd.read_csv(dir_work + 'OGindelStat.lst', sep = '\t', header = 0)
df_indel['OGid'] = df_indel['OGid'].apply(str)
df_indel['indel'] = np.where(df_indel['species'] == 'pma_1', 'del', 'ins')
#df_indel['OGid_no'] = df_indel.groupby(['OGid', 'indel'])['start'].rank('max').astype(int).astype(str)
df_indel['nw_OGid'] = df_indel['OGid'] + '.' + df_indel['indel'] + '.' + df_indel['start'].apply(str) + '.' + df_indel['end'].apply(str)
indel_mid = df_indel[df_indel['pos'] == 'mid']['nw_OGid'].tolist()


# make bed files
lst_sp = ['pma_1', 'tgu_1', 'fal_1', 'cco_1']
buffer = 300

for index, row in df_indel.iterrows():
    id_no = row['nw_OGid']
    if id_no in indel_mid:
        indel_up = row['start'] - buffer
        indel_dwn = row['end'] + buffer
        if indel_up < 0: 
            indel_up = 0
        if indel_dwn > row['OGlen']: 
            indel_dwn = row['OGlen']
        for sp in lst_sp:
            up = sp + '\t' + str(indel_up) + '\t' + str(row['start']) + '\n'
            dwn = sp + '\t' + str(row['end'])  + '\t' + str(indel_dwn) + '\n'
            with open(dir_work + 'bed_surr/' + id_no + '.up.bed', 'a+') as outUp, open(dir_work + 'bed_surr/' + id_no + '.dwn.bed', 'a+') as outdwn: # surr:surrounding
                outUp.write(up)
                outdwn.write(dwn)
    else:
        if row['start'] == 0:
            indel_dwn = row['end'] + buffer
            if indel_dwn > row['OGlen']: 
                indel_dwn = row['OGlen']
            for sp in lst_sp:
                dwn = sp + '\t' + str(row['end'])  + '\t' + str(indel_dwn) + '\n'
                with open(dir_work + 'bed_surr/' + id_no + '.dwn.bed', 'a+') as outdwn:
                    outdwn.write(dwn)
        else:
            for sp in lst_sp:
                indel_up = row['start'] - buffer
                up = sp + '\t' + str(indel_up) + '\t' + str(row['start']) + '\n'
                with open(dir_work + 'bed_surr/' + id_no + '.up.bed', 'a+') as outUp: # surr:surrounding
                    outUp.write(up)

# get alignment nt files
lst_bed = glob.glob(dir_work + 'bed_surr/*.bed')
os.system('mkdir -p ' + dir_work + 'cds/')
os.system('rm ' + dir_work + 'seq/*/*.fai')

for bed in lst_bed:
    header = re.sub('^.*/', '', bed)
    #print(header)
    tag = header.split('.')
    #pos = tag[2]
    id_no = tag[0]
    id_in = re.sub(r'\.bed', '', header)
    print(id_in)
    os.system('bedtools getfasta -fi ' + dir_work + 'seq/' + id_no + '/' + id_no + '_nostp.fasta -bed ' + dir_work + 'bed_surr/' + id_in + '.bed > ' + dir_work + 'cds/' + id_in + '.fasta')

#os.system('sed -i \'s/:[0-9]*-[0-9]*//g\'' + dir_work + ' cds/*.fasta')


# Make dirs for different dataset
subdirs = ['all/', 'ins/', 'del/']  # List of subdirectories to create
for subdir in subdirs:
    os.makedirs(os.path.join(dir_out, subdir), exist_ok=True)

# import data
lst_cds_all = glob.glob(dir_work + r'cds/*.fasta')
lst_cds_ins = glob.glob(dir_work + r'cds/*ins*.fasta')
lst_cds_del = glob.glob(dir_work + r'cds/*del*.fasta')

def aln_from_file(filepath):
    alignment_dict = defaultdict(str)
    # Read the alignment file (can handle FASTA, CLUSTAL, etc.)
    alignment = AlignIO.read(filepath, "fasta")  # Adjust format if needed
    # Store sequences in the dictionary
    for record in alignment:
        alignment_dict[record.id] = str(record.seq)
    return alignment_dict

def import_aln(filepaths):
    all_alignments = {}
    for filepath in filepaths:
        fname = os.path.basename(filepath)[:-6]
        alignment = aln_from_file(filepath)
        all_alignments[fname] = alignment
    return all_alignments

# Set dynamic lengths based on the longest sequence found
def process_aln(lst_aln, buffer, seq_types=['up', 'dwn']):
    codon_concat = defaultdict(lambda: defaultdict(str))
    # Create keys for each seq_type (up, dwn)
    for seq_type in seq_types:
        for i in range(1, int(buffer / 3) + 1):
            key = f'{seq_type}{i}'  # Generate keys like up1, up2, ..., dwn100
    # Iterate over the sequences from lst_aln
    for seq_label, sequences in lst_aln.items():
        group = defaultdict(str)  # Group sequences by species
        for seq_type in seq_types:
            if seq_label.endswith(f".{seq_type}"):
                # Concatenate sequences for each species
                for spe, seq in sequences.items():
                    speID = spe.split(':')[0]  # Extract species ID (e.g., pma_1)
                    group[speID] += seq  # Concatenate sequences for each species
                # Process concatenated sequences
                for species, concat_seq in group.items():
                    # Split the concatenated sequence into codons (chunks of 3)
                    codons = [concat_seq[j:j+3] for j in range(0, len(concat_seq), 3)]
                    # Loop through the 100 keys and populate codon_concat
                    for i in range(1, int(buffer / 3) + 1):
                        key = f'{seq_type}{i}'  # Define key (up1, up2, ..., dwn100)
                        if i <= len(codons):  # Ensure there are enough codons for the key
                            # Append the codon to the appropriate species in codon_concat
                            codon_concat[key][species] += codons[i - 1]
    return codon_concat

# output to file
def codon2fas(codon_concat, dir_out):
    # Loop through the codon_concat dictionary
    for codon_name, sequences in codon_concat.items():
        # Open a new file in the directory 'all/' with the codon_name as the filename
        file_path = f"{dir_out}/{codon_name}.fasta"
        with open(file_path, "w") as f:
            # Loop through each species sequence in sequences
            for seq_name, sequence in sequences.items():
                # Write the header in FASTA format
                f.write(f">{seq_name}\n")
                # Write the sequence in blocks of 60 characters per line (typical FASTA format)
                for i in range(0, len(sequence), 60):
                    f.write(f"{sequence[i:i+60]}\n")

aln = [(lst_cds_all, dir_out + 'all/'), (lst_cds_ins, dir_out + 'ins/'), (lst_cds_del, dir_out + 'del/')]
for file_list, output_dir in aln:
    # Import alignments
    lst_aln = import_aln(file_list)
    # Process alignments
    codon_concat = process_aln(lst_aln, buffer)
    # Output to FASTA files
    codon2fas(codon_concat, output_dir)

"""
# Fasta files to phylip
def fas2phylip(input_fasta, output_phylip):
    # Read the input FASTA file and convert to PHYLIP format
    with open(input_fasta, "r") as fasta_file:
        alignments = AlignIO.read(fasta_file, "fasta")
    # Write the alignments in PHYLIP format
    with open(output_phylip, "w") as phylip_file:
        AlignIO.write(alignments, phylip_file, "phylip")
"""

lst_fa = glob.glob(dir_out + '*/*.fasta')
for fa in lst_fa:
    #print(fa)
    path = re.sub(r'\.fasta', '', fa)
    #print(path)
    os.system('perl /work/myucchen/software/Fasta2Phylip.pl ' + fa + ' ' + path + '.phylip')

os.system('sed -i \'s/\t/  /g\' ' + dir_out + '*/*.phylip')

# paml ctl
base_in_dir = '/work/myucchen/indel/indel.03.02/inputs/{}/'
base_out_dir = '/work/myucchen/indel/indel.03.02/{}/m{}/out/'
base_sh_dir = '/work/myucchen/indel/indel.03.02/sh/{}/m{}/'
dat_path = '/work/myucchen/indel/indel.02.07/dat.path.list'
tree_path = '/work/myucchen/indel/indel.02.07/species.tre'
perl_script = '/work/myucchen/script/cmd_codeml.3.pl'
codeml_exe = '/work/myucchen/software/paml-4.10.7/bin/codeml'

# Define the groups and models
groups = ['del', 'ins', 'all']
models = {
    'm0': '0',
    'm1': '1',
    'm2': '2'
}


# Loop through each group and model
for group in groups:
    for model_name, model_value in models.items():
        # Construct the directories for input, output, and shell scripts
        in_dir = base_in_dir.format(group)
        out_dir = base_out_dir.format(group, model_name[1])  # 'm0' -> '0', 'm1' -> '1', etc.
        sh_dir = base_sh_dir.format(group, model_name[1])
        # Construct the Perl command
        perl_cmd = (
            f'perl {perl_script} --exe={codeml_exe} '
            f'--in_dir={in_dir} --in_file_ext=phylip '
            f'--out_dir={out_dir} --dat={dat_path} --tree={tree_path} '
            f'--sh_dir={sh_dir} --seqtype="1" --CodonFreq="2" '
            f'--model="{model_value}" --NSsites="0" --icode="0" '
            f'--n_job="10" --debug="1" --omega="2" --fix_omega="0"'
        )
        # Print the command for each combination of group and model
        print(f"Running command for group {group}, model {model_name}:")
        print(perl_cmd)
        os.system(perl_cmd)
        print()


# paml execute
cmd_no = 10 #int(buffer / 3)
dir_base = '/work/myucchen/indel/indel.03.02/sh/{}/{}/'
groups = ['del', 'ins', 'all']
models = ['m0', 'm1', 'm2']

for group in groups:
    for model in models:
        dir_sh = dir_base.format(group, model)
        os.system('mkdir -p ' + dir_sh)
        for i in range(1, cmd_no + 1):
            with open(dir_sh + 'exe{}.sh'.format(i), 'w') as exe_out:
                exe_out.write('#!/bin/bash\n')
                # sbatch setting
                set_sbatch = '#SBATCH --job-name paml{}\n#SBATCH --partition=med\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'exe{}.out\n#SBATCH -e '.format(i) + dir_sh + 'exe{}.err\n\n'.format(i)
                exe_out.write(set_sbatch)
                run = 'srun -n 1 -c 1 ' + dir_sh + 'job{}.sh'.format(i) + ' &\n\n'
                exe_out.write(run)
                exe_out.write('wait\n')
        with open(dir_sh + 'paml.sh', 'w') as exe_out:
            exe_out.write('#!/bin/bash\n')
            for i in range(1, cmd_no + 1):
                exe_out.write('sbatch ' + dir_sh + 'exe{}.sh'.format(i) + ' &\n')
            exe_out.write('wait\n')
        # Change directory permissions and run the paml.sh script
        os.chdir(dir_sh)
        os.system('chmod -R 755 ' + dir_sh[:-11])
        os.system('nohup bash ' + dir_sh + 'paml.sh &')


#####CONTROL REGION
# paml for control region
dir_work = '/work/myucchen/indel/indel.03.01/'
dir_seq = dir_work + 'seq/'
dir_in = dir_work + 'bed_ctrl/'
dir_out = '/work/myucchen/indel/indel.03.02/ctrl/'

os.system('mkdir -p ' + dir_out + 'cds/')
os.system('mkdir -p ' + dir_out + 'inputs/')
os.system('rm ' + dir_work + 'seq/*/*.fai')

lst_bed = glob.glob(dir_in + '*.bed')

for bed in lst_bed:
    header = re.sub('^.*/', '', bed)
    #print(header)
    tag = header.split('.')
    #print(tag)
    id_no = tag[0]
    os.system('bedtools getfasta -fi ' + dir_work + 'seq/' + id_no + '/' + id_no + '_nostp.fasta -bed ' + bed + ' > ' + dir_out + 'cds/' + id_no + '.fasta')

# Make dirs for different dataset
subdirs = ['all/']  # List of subdirectories to create
for subdir in subdirs:
    os.makedirs(os.path.join(dir_out, 'inputs', subdir), exist_ok=True)

# import data
lst_cds_all = glob.glob(dir_out + r'cds/*.fasta')
aln = [(lst_cds_all, dir_out + 'inputs/all/')]

def process_aln_ctrl(lst_aln, buffer):
    codon_concat = defaultdict(lambda: defaultdict(str))
    for i in range(1, int(buffer / 3) + 1):
        key = f'pos{i}'  # Generate keys like pos1, pos2, ..., pos100
    # Iterate over the sequences from lst_aln
    for seq_label, sequences in lst_aln.items():
        group = defaultdict(str)  # Group sequences by species
        # Concatenate sequences for each species
        for spe, seq in sequences.items():
            speID = spe.split(':')[0]  # Extract species ID (e.g., pma_1)
            group[speID] += seq  # Concatenate sequences for each species
        # Process concatenated sequences
        for species, concat_seq in group.items():
            # Split the concatenated sequence into codons (chunks of 3)
            codons = [concat_seq[j:j + 3] for j in range(0, len(concat_seq), 3)]
            # Loop through the 100 keys and populate codon_concat
            for i in range(1, int(buffer / 3) + 1):
                key = f'pos{i}'  # Define key (pos1, pos2, ..., pos100)
                if i <= len(codons):  # Ensure there are enough codons for the key
                    # Append the codon to the appropriate species in codon_concat
                    codon_concat[key][species] += codons[i - 1]
    return codon_concat

for file_list, output_dir in aln:
    # Import alignments
    lst_aln = import_aln(file_list)
    # Process alignments
    codon_concat = process_aln_ctrl(lst_aln, buffer)
    # Output to FASTA files
    codon2fas(codon_concat, output_dir)


lst_fa = glob.glob(dir_out + 'inputs/*/*.fasta')
for fa in lst_fa:
    #print(fa)
    path = re.sub(r'\.fasta', '', fa)
    #print(path)
    os.system('perl /work/myucchen/software/Fasta2Phylip.pl ' + fa + ' ' + path + '.phylip')

os.system('sed -i \'s/\t/  /g\' ' + dir_out + '*/*/*.phylip')

# paml ctl
base_in_dir = '/work/myucchen/indel/indel.03.02/ctrl/inputs/{}/'
base_out_dir = '/work/myucchen/indel/indel.03.02/ctrl/{}/m{}/out/'
base_sh_dir = '/work/myucchen/indel/indel.03.02/ctrl/sh/{}/m{}/'
dat_path = '/work/myucchen/indel/indel.02.07/dat.path.list'
tree_path = '/work/myucchen/indel/indel.02.07/species.tre'
perl_script = '/work/myucchen/script/cmd_codeml.3.pl'
codeml_exe = '/work/myucchen/software/paml-4.10.7/bin/codeml'

# Define the groups and models
groups = ['all']
models = {
    'm0': '0',
    'm1': '1',
    'm2': '2'
}

# Loop through each group and model
for group in groups:
    for model_name, model_value in models.items():
        # Construct the directories for input, output, and shell scripts
        in_dir = base_in_dir.format(group)
        out_dir = base_out_dir.format(group, model_name[1])  # 'm0' -> '0', 'm1' -> '1', etc.
        sh_dir = base_sh_dir.format(group, model_name[1])
        # Construct the Perl command
        perl_cmd = (
            f'perl {perl_script} --exe={codeml_exe} '
            f'--in_dir={in_dir} --in_file_ext=phylip '
            f'--out_dir={out_dir} --dat={dat_path} --tree={tree_path} '
            f'--sh_dir={sh_dir} --seqtype="1" --CodonFreq="2" '
            f'--model="{model_value}" --NSsites="0" --icode="0" '
            f'--n_job="100" --debug="1" --omega="2" --fix_omega="0"'
        )
        # Print the command for each combination of group and model
        print(f"Running command for group {group}, model {model_name}:")
        print(perl_cmd)
        os.system(perl_cmd)
        print()


# paml execute
cmd_no = int(buffer / 3)
dir_base = '/work/myucchen/indel/indel.03.02/ctrl/sh/{}/{}/'
groups = ['all']
models = ['m0', 'm1', 'm2']

for group in groups:
    for model in models:
        dir_sh = dir_base.format(group, model)
        os.system('mkdir -p ' + dir_sh)
        for i in range(1, cmd_no + 1):
            with open(dir_sh + 'exe{}.sh'.format(i), 'w') as exe_out:
                exe_out.write('#!/bin/bash\n')
                # sbatch setting
                set_sbatch = '#SBATCH --job-name paml{}\n#SBATCH --partition=med\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'exe{}.out\n#SBATCH -e '.format(i) + dir_sh + 'exe{}.err\n\n'.format(i)
                exe_out.write(set_sbatch)
                run = 'srun -n 1 -c 1 ' + dir_sh + 'job{}.sh'.format(i) + ' &\n\n'
                exe_out.write(run)
                exe_out.write('wait\n')
        with open(dir_sh + 'paml.sh', 'w') as exe_out:
            exe_out.write('#!/bin/bash\n')
            for i in range(1, cmd_no + 1):
                exe_out.write('sbatch ' + dir_sh + 'exe{}.sh'.format(i) + ' &\n')
            exe_out.write('wait\n')
        # Change directory permissions and run the paml.sh script
        os.chdir(dir_sh)
        os.system('chmod -R 755 ' + dir_sh[:-11])
        os.system('nohup bash ' + dir_sh + 'paml.sh &')



















lst_fa = glob.glob(dir_out + '*.fasta')
for fa in lst_fa:
    #print(fa)
    path = re.sub(r'\.fasta', '', fa)
    #print(path)
    os.system('perl /work/myucchen/software/Fasta2Phylip.pl ' + fa + ' ' + path + '.phylip')

os.system('sed -i \'s/\t/  /g\' ' + dir_out + '*.phylip')



# paml ctl
base_in_dir = '/work/myucchen/indel/indel.03.02/ctrl/inputs/{}/'
base_out_dir = '/work/myucchen/indel/indel.03.02/ctrl/{}/m{}/out/'
base_sh_dir = '/work/myucchen/indel/indel.03.02/ctrl/sh/{}/m{}/'
dat_path = '/work/myucchen/indel/indel.02.07/dat.path.list'
tree_path = '/work/myucchen/indel/indel.02.07/species.tre'
perl_script = '/work/myucchen/script/cmd_codeml.3.pl'
codeml_exe = '/work/myucchen/software/paml-4.10.7/bin/codeml'

# Define the groups and models
groups = ['all']
models = {
    'm0': '0',
    'm1': '1',
    'm2': '2'
}

# Loop through each group and model
for group in groups:
    for model_name, model_value in models.items():
        # Construct the directories for input, output, and shell scripts
        in_dir = base_in_dir.format(group)
        out_dir = base_out_dir.format(group, model_name[1])  # 'm0' -> '0', 'm1' -> '1', etc.
        sh_dir = base_sh_dir.format(group, model_name[1])
        # Construct the Perl command
        perl_cmd = (
            f'perl {perl_script} --exe={codeml_exe} '
            f'--in_dir={in_dir} --in_file_ext=phylip '
            f'--out_dir={out_dir} --dat={dat_path} --tree={tree_path} '
            f'--sh_dir={sh_dir} --seqtype="1" --CodonFreq="2" '
            f'--model="{model_value}" --NSsites="0" --icode="0" '
            f'--n_job="100" --debug="1" --omega="2" --fix_omega="0"'
        )
        # Print the command for each combination of group and model
        print(f"Running command for group {group}, model {model_name}:")
        print(perl_cmd)
        os.system(perl_cmd)
        print()


# paml execute
cmd_no = int(buffer / 3)
dir_base = '/work/myucchen/indel/indel.03.02/sh/{}/{}/'
groups = ['del', 'ins', 'all']
models = ['m0', 'm1', 'm2']

for group in groups:
    for model in models:
        dir_sh = dir_base.format(group, model)
        os.system('mkdir -p ' + dir_sh)
        for i in range(1, cmd_no + 1):
            with open(dir_sh + 'exe{}.sh'.format(i), 'w') as exe_out:
                exe_out.write('#!/bin/bash\n')
                # sbatch setting
                set_sbatch = '#SBATCH --job-name paml{}\n#SBATCH --partition=med\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'exe{}.out\n#SBATCH -e '.format(i) + dir_sh + 'exe{}.err\n\n'.format(i)
                exe_out.write(set_sbatch)
                run = 'srun -n 1 -c 1 ' + dir_sh + 'job{}.sh'.format(i) + ' &\n\n'
                exe_out.write(run)
                exe_out.write('wait\n')
        with open(dir_sh + 'paml.sh', 'w') as exe_out:
            exe_out.write('#!/bin/bash\n')
            for i in range(1, cmd_no + 1):
                exe_out.write('sbatch ' + dir_sh + 'exe{}.sh'.format(i) + ' &\n')
            exe_out.write('wait\n')
        # Change directory permissions and run the paml.sh script
        os.chdir(dir_sh)
        os.system('chmod -R 755 ' + dir_sh[:-11])
        os.system('nohup bash ' + dir_sh + 'paml.sh &')



"""
# Build the dict for the output
codon_concat = defaultdict(lambda: defaultdict(str))
for seq_type in ['up', 'dwn']:
    # For each seq_type (up/dwn), create 100 keys (up1 to up100, dwn1 to dwn100)
    for i in range(1, int(buffer/3+1)):
        key = f'{seq_type}{i}'  # Generate key like up1, up2, ..., dwn100
        # For each species, concatenate their sequences across regions
        key = f'{seq_type}{i}'

for seq_label, sequences in lst_aln.items():
    group = defaultdict(str)  # Group sequences by species
    for seq_type in ['up', 'dwn']:
        if seq_label.endswith(f".{seq_type}"):
            # Go through each species in the sequences
            for spe, seq in sequences.items():
                speID = spe.split(':')[0]  # Extract species ID (e.g., pma_1)
                group[speID] += seq  # Concatenate sequences for each species
            for species, concat_seq in group.items():
                #print(group)  # Debug: print grouped sequences
                # Split the concatenated sequence into chunks of 3 (codons)
                codons = [concat_seq[j:j+3] for j in range(0, len(concat_seq), 3)]
                # Check if we have enough codons to fill the current key (up to 100)
                for i in range(1, int(buffer/3+1)):  # Ensure we loop through the 100 keys again
                    key = f'{seq_type}{i}'  # Define the key for codon_concat
                    if i <= len(codons):  # If we have enough codons for this key
                        # Append the codon to the appropriate species in codon_concat
                        codon_concat[key][species] += codons[i-1]  # Append the codon

# output to file
for codon_name, sequences in codon_concat.items():
    # Open a new file with the codon_name as filename
    with open(dir_out + 'all/' + f"{codon_name}.fasta", "w") as f:
        for seq_name, sequence in sequences.items():
            # Write the header and sequence in FASTA format
            f.write(f">{seq_name}\n")
            # Write the sequence in blocks of 60 characters per line (typical FASTA format)
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")



cmd_no = int(buffer / 3)
dir_sh = '/work/myucchen/indel/indel.02.09/sh/del/m0/'

for i in range(1, cmd_no + 1):
    with open(dir_sh + 'exe{}.sh'.format(i), 'w') as exe_out:
        exe_out.write('#!/bin/bash\n')
        # sbatch setting
        set_sbatch = '#SBATCH --job-name paml{}\n#SBATCH --partition=med\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'exe{}.out\n#SBATCH -e '.format(i) + dir_sh + 'exe{}.err\n\n'.format(i)
        exe_out.write(set_sbatch)
        run = 'srun -n 1 -c 1 ' + dir_sh + 'job{}.sh'.format(i) + ' &\n\n'
        exe_out.write(run)
        exe_out.write('wait\n')

with open(dir_sh + 'paml.sh', 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(1, cmd_no + 1):
        exe_out.write('sbatch ' + dir_sh + 'exe{}.sh'.format(i) + ' &\n')
    exe_out.write('wait\n')

os.chdir(dir_sh)
os.system('chmod -R 755 ' + dir_sh[:-11])
os.system('nohup bash ' + dir_sh + 'paml.sh &')

"""
