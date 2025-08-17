### alignment and ancestral reconstruction
### indel_aln_recon.2.py
### v2 2024/07/11
### utilize list from orthofinder
### v1 2024/06/07
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
from Bio.Seq import Seq
from Bio import AlignIO

dir_work = '/work/myucchen/indel/indel.02.04/'
dir_in = '/work/myucchen/indel/indel.02.03/aa/'
dir_out = dir_work + 'seq/'
dir_sh = dir_work + 'sh/'
indelMp = '/work/myucchen/software/indelMaP/src/indelMaP_ASR.py'
tre = dir_in[:-3] + 'species.tre'
lst_spe = '/work/myucchen/indel/indel.02.01/species.lst'
lst_sngl_cp = '/work/myucchen/indel/indel.02.02/orthofinder/Results_Aug16/Orthogroups/Orthogroups_SingleCopyOrthologues.txt'
fa_ext = 'fasta'
cmd_no = 100

# split list into sublist
def chunkify(total_item, item_in_lst):
    return [total_item[i::item_in_lst] for i in list(range(item_in_lst))]

os.system('mkdir -p ' + dir_out)
os.system('mkdir -p ' + dir_sh)

# alignment
lst_fa = glob.glob(dir_in + '*.' + fa_ext)
# single-copy orthogroups
with open(lst_sngl_cp, 'r+') as fileIn:
    sngl_cp = [line.strip() for line in fileIn]

sngl_cp = [re.sub(r'OG0+(.+)', r'\1', og) for og in sngl_cp]
lst_aln = list(chunkify(sngl_cp, cmd_no))

"""
for fa in lst_aln:
    file_name = re.sub('^.*/', '', fa)
    file_name = re.sub('.fasta', '', file_name)
    os.system('prank -d=' + fa + ' -t=' + tre + ' -showall +F -o=' + dir_out + file_name)
"""

# alignment
os.system('mkdir -p ' + dir_sh + 'aln')
i = 0
for index, line in enumerate(lst_aln):
    with open(dir_sh + 'aln/job_aln{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/home/myucchen/miniconda3/bin:$PATH\nsource activate base\n')
        for og_id in lst_aln[i]:
            #file_name = re.sub('^.*/', '', og_id)
            #file_name = re.sub('.fasta', '', file_name)
            #print(file_name)
            print(og_id)
            mkdir = 'mkdir -p ' + dir_out + og_id + '/' + '\ncd ' + dir_out + og_id + '\n'
            aln = 'prank -d=' + dir_in + og_id + '.fasta -t=' + tre + ' -showall +F -o=' + dir_out + og_id + '/' + og_id + '\n'
            #asr_cmd = 'python ' + asr + ' --msa_file ' +  dir_out + file_name + '/' + file_name + '.best.fas --tree_file ' + tre + '  --alphabet Protein --output_file ' +  dir_out + file_name + '/' + file_name + '\n'
            bash_out.write(mkdir + aln)
        i += 1

for i in range(0, cmd_no):
    with open(dir_sh + 'aln/cmd_aln{}.sh'.format(i), 'w') as exe_out:
        exe_out.write('#!/bin/bash\n')
        # sbatch setting
        set_sbatch = '#SBATCH --job-name aln{}\n#SBATCH --partition=long\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'aln/cmd_aln{}.out\n#SBATCH -e '.format(i) + dir_sh + 'aln/cmd_aln{}.err\n\n'.format(i)
        exe_out.write(set_sbatch)
        cal = 'srun -n 1 -c 1 ' + dir_sh + 'aln/job_aln{}.sh'.format(i) + ' &\n\n'
        exe_out.write(cal)
        exe_out.write('wait\n')

with open(dir_sh + 'aln/exe_aln.sh', 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('sbatch ' + dir_sh + 'aln/cmd_aln{}.sh'.format(i) + ' &\n')
    exe_out.write('wait\n')

os.chdir(dir_sh)
os.system('chmod -R 777 ' + dir_work)
os.system('nohup bash ' + dir_sh + 'aln/exe_aln.sh &')
# 12421 genes

# find indel potentially not in frame
# check the cds sequences can not be divided by 3
not_inframe = []
excl_og = []
for seqID in sngl_cp:
    with open(dir_in[:-3] + 'cds/' + seqID + '.fasta', 'r') as inSeq:
        for seq in SeqIO.parse(inSeq, 'fasta'):
            #print(seq.name)
            #print(len(seq.seq))
            if len(seq.seq)%3 != 0:
                #print(seqID + '\n' + seq.name)
                not_inframe.append(seqID + ':' + seq.name + ':' + str(len(seq.seq)))
                excl_og.append(seqID)

len(excl_og)
# 33
excl_og = list(set(excl_og))
len(excl_og)
# 27

with open(dir_work + 'not_inframe.lst', 'w+') as outLst:
    for item in not_inframe:
        outLst.write(item + '\n')

# exclude orthogroups with sequences (in any species) not in frame
lst_asr = [x for x in sngl_cp if x not in excl_og]
lst_asr =  list(chunkify(lst_asr, cmd_no))

# make codon alignment and ancestral reconstruction
os.system('mkdir -p ' + dir_sh + 'asr')
i = 0
for index, line in enumerate(lst_asr):
    with open(dir_sh + 'asr/job_asr{}.sh'.format(index), 'w') as bash_out:
        bash_out.write('#!/bin/bash\nexport PATH=/home/myucchen/miniconda3/bin:$PATH\nsource activate base\n')
        for og_id in lst_asr[i]:
            print(og_id)
            to_dir = 'cd ' + dir_out + og_id + '\n'
            srt_cds = 'seqkit sort -n ' + dir_in[:-3] + 'cds/' + og_id + '.fasta -o ' + og_id + '.cds.sort.fas\n'
            srt_aa = 'seqkit sort -n ' + og_id + '.best.fas -o ' + og_id + '.aa.sort.fas\n'
            codon = 'pal2nal.pl ' + og_id + '.aa.sort.fas ' + og_id + '.cds.sort.fas -output fasta -codontable 1 > ' + og_id + '.codon\n' 
            asr = 'python ' + indelMp + ' --msa_file ' +  dir_out + og_id + '/' + og_id + '.codon --tree_file ' + tre + '  --alphabet DNA --output_file ' +  dir_out + og_id + '/' + og_id + '\n'
            rm = 'rm ' + og_id + '.aa.sort.fas ' + og_id + '.cds.sort.fas\n'
            bash_out.write(to_dir + srt_cds + srt_aa + codon + asr + rm)
        i += 1


for i in range(0, cmd_no):
    with open(dir_sh + 'asr/cmd_asr{}.sh'.format(i), 'w') as exe_out:
        exe_out.write('#!/bin/bash\n')
        # sbatch setting
        set_sbatch = '#SBATCH --job-name asr{}\n#SBATCH --partition=long\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem=8G\n#SBATCH -o '.format(i) + dir_sh + 'asr/cmd_asr{}.out\n#SBATCH -e '.format(i) + dir_sh + 'asr/cmd_asr{}.err\n\n'.format(i)
        exe_out.write(set_sbatch)
        cal = 'srun -n 1 -c 1 ' + dir_sh + 'asr/job_asr{}.sh'.format(i) + ' &\n\n'
        exe_out.write(cal)
        exe_out.write('wait\n')

with open(dir_sh + 'asr/exe_asr.sh', 'w') as exe_out:
    exe_out.write('#!/bin/bash\n')
    for i in range(0, cmd_no):
        exe_out.write('sbatch ' + dir_sh + 'asr/cmd_asr{}.sh'.format(i) + ' &\n')
    exe_out.write('wait\n')

os.chdir(dir_sh + 'asr/')
os.system('chmod -R 777 ' + dir_work)
os.system('nohup bash ' + dir_sh + 'asr/exe_asr.sh &')



# find the stop codon in the files
# not in frame list
df_notFrame = pd.read_csv(dir_work + 'not_inframe.lst', sep = ':', names = ['OG', 'seqID', 'seqLen'])
lst_notFrame = df_notFrame['OG'].astype(str).tolist()
lst_notFrame = list(set(lst_notFrame))

codon_stop = ["TAG", "TGA", "TAA"]
lst_og = glob.glob(dir_work + 'seq/*')

for dir_og in lst_og:
    os.chdir(dir_og)
    id_no = re.sub('^.*/', '', dir_og)
    if id_no in lst_notFrame: continue
    print('ID: ' + id_no)
    os.system('cat ' + id_no + '.codon ' + id_no + '_internal_ancestral_reconstruction.fas > ' +  id_no + '.' + fa_ext)


df_instp = pd.DataFrame() # internal stop codon
df_stp = pd.DataFrame() # all stop codon

lst_fa = glob.glob(dir_work + 'seq/*/*.fasta')

for fa in lst_fa:
    og_id = re.sub('^.*/', '', fa)
    og_id = re.sub('.fasta','',og_id)
    print(og_id)
    for seq in SeqIO.parse(fa, "fasta"):
        tmp = list(seq.seq)
        #print(seq.seq)
        # Iterate over the sequence in steps of 3 to identify codons
        for index in range(0, len(seq.seq), 3):
            codon = str(seq.seq[index:index+3])  # Convert the Seq object to string
            if codon in codon_stop:
                if index+3 != len(seq.seq):
                    stp = pd.DataFrame({'og_id': [og_id], 'header': [seq.name], 'start': [index], 'end': [index+3]})
                    df_instp = pd.concat([df_instp, stp])
                else:
                    stp = pd.DataFrame({'og_id': [og_id], 'header': [seq.name], 'start': [index], 'end': [index+3]})
                    df_stp = pd.concat([df_stp, stp])


#df = df.drop('column_name', axis=1)
df_instp = df_instp.reset_index(drop = True)
df_instp.to_csv(dir_work + 'internal_stop_codon.lst', sep = '\t', index = None)
df_instp_pos = df_instp[['og_id', 'start', ]]
df_instp_pos = df_instp_pos.drop_duplicates(['og_id', 'start'], keep = 'first').reset_index(drop = True) # 1336 internal stop codon
#lst_instp_pos = df_instp_pos[['og_id']].drop_duplicates(['og_id'], keep = 'first') # 939 ortholog groups

df_stp = df_stp.reset_index(drop = True)
df_stp_pos = df_stp[['og_id', 'start']] 
df_stp_pos = df_stp_pos.drop_duplicates(['og_id', 'start'], keep = 'first').reset_index(drop = True) # 12382 ending stop codon
# df_stp_pos[['og_id']].drop_duplicates(['og_id'], keep = 'first') # 12382 rows
# 12 orthogroups don't have ending stop codon

df_stp_all = pd.concat([df_instp_pos, df_stp_pos]) # 13718 stop codon in total


lst_fa = glob.glob(dir_work + 'seq/*/*.fasta')

for fa in lst_fa:
    og_id = re.sub('^.*/', '', fa)
    og_id = re.sub('.fasta','',og_id)
    print(og_id)
    all_stp = df_stp_all.loc[df_stp_all['og_id'] == og_id,:].reset_index(drop = True)
    start = all_stp['start'].tolist()
    #print(start)
    seqs = list(SeqIO.parse(fa, 'fasta'))
    seqs_new = []
    for seq in seqs:
        tmp = list(seq.seq)
        start_sorted = sorted(start, reverse=True)
        for index in start_sorted:
            #print(tmp[index:index+3])
            tmp[index:index+3] = [] # remove the codon from temporary sequences
        #print((''.join(tmp)))
        seq.seq = Seq(''.join(tmp))
        seqs_new.append(seq)
    SeqIO.write(seqs_new, dir_work + 'seq/' + og_id + '/' + og_id + '_nostp.fasta', 'fasta')
