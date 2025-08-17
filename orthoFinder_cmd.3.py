### Prepare seqs for Orthofinder
### orthoFinder_cmd.3.py
### v3 2025/03/14
### For ensembl dataset, longest isoform dataset
### v2 2025/03/14
### For ensembl dataset
### v1 2024/05/08
### Prepare files for orthofinder
### @author:Chen Yu Chi

import os,sys
import glob
import re
import csv
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO

dir_Seqin = '/vol/storage/indel/indel.01.04/cds/' # organize files for aa
dir_Seqout = '/vol/storage/indel/indel.01.04/cds_ortho/' # organize files for aa
file_ext = 'fasta'
file_species = '/vol/storage/indel/indel.01.04/species.list'
file_tre = '/vol/storage/indel/indel.01.01/tre/363-avian-2020-phast_prn.tre'
dir_out = 'orthofinder'
#lst_species = ['PARMAJ', 'TAEGUT', 'FICALB']

os.system('mkdir -p ' + dir_Seqout)

# change file name
lst_fa = glob.glob(dir_Seqin + '*.' + file_ext)
df_species = pd.read_csv(file_species, sep = '\t', names = ['species', 'id', 'abbr'])
for fa in lst_fa:
    last_slash_index = fa.rfind('/')
    id_file = fa[last_slash_index + 1:]
    # id_file = id_file[:15] # ncbi
    id_file = re.sub(r'\..*.fasta', '', id_file) # ensembl
    # print(id_file)
    id_sp = id_file
    id_abbr = df_species[df_species['species'] == id_file]['abbr'].values[0]
    print(id_abbr)
    os.system('cp ' + fa + ' ' + dir_Seqout + id_sp + '.fasta')


# parsing ensembl fasta header information to dataframe
"""
def parse_data(entry):
    elements = entry.split()
    parsed = {
        "Transcript_ID": elements[0],
        "Type": elements[1],
        "Assembly": elements[2].split(":")[1],
        "Chromosome": elements[2].split(":")[2],
        "Start": elements[2].split(":")[3],
        "End": elements[2].split(":")[4],
        "Strand": elements[2].split(":")[5],
        "Gene_ID": elements[3].split(":")[1],
        "Gene_Biotype": elements[4].split(":")[1],
        "Transcript_Biotype": elements[5].split(":")[1],
        "Gene_Symbol": None,
        "Description": None
    }
    for element in elements[6:]:
        if element.startswith("gene_symbol:"):
            parsed["Gene_Symbol"] = element.split(":")[1]
        elif element.startswith("description:"):
            parsed["Description"] = " ".join(elements[elements.index(element):]).replace("description:", "")
            break
    return parsed
"""

def parse_data(entry):
    elements = entry.split()
    parsed = {
        "Transcript_ID": None,
        "Gene_ID": None,
        "Seq_ID": None,
        "Type": None
    }
    for element in elements:
        if element.startswith("transcript:"):
            parsed["Transcript_ID"] = element.split(":")[1]
        elif element.startswith("gene=gene:"):
            parsed["Gene_ID"] = element.split(":")[1]
        elif element.startswith("seq_id="):
            parsed["Seq_ID"] = element.split("=")[1]
        elif element.startswith("type="):
            parsed["Type"] = element.split("=")[1]
    return parsed


lst_faInfo = []
lst_fa_h = glob.glob(dir_Seqout + '*.fasta')
for fa in lst_fa_h:
    last_slash_index = fa.rfind('/')
    id_file = fa[last_slash_index + 1:]
    id_file = re.sub(r'\.fasta', '', id_file)
    # print(id_file)
    with open(fa, 'r+') as faIn:
        for record in SeqIO.parse(faIn, "fasta"):
            print(record.description)
            parsed_entry = parse_data(record.description)
            parsed_entry["species"] = id_file  # Add species column
            lst_faInfo.append(parsed_entry)

df_faInfo = pd.DataFrame(lst_faInfo)
#df_faInfo.to_csv(re.sub(r'^(.*\/)[^\/]+\/?$', r'\1', dir_Seqout) + 'headerInfo.lst', sep = '\t', na_rep= "NA", index = None, header = ['transcriptID', 'type', 'assembly', 'chr', 'start', 'end', 'strand', 'geneID', 'gene_biotyp', 'transcript_biotyp', 'gen_sym', 'description', 'species'])
df_faInfo.to_csv(re.sub(r'^(.*\/)[^\/]+\/?$', r'\1', dir_Seqout) + 'headerInfo.lst', sep = '\t', na_rep= "NA", index = None, header = ['transcriptID', 'Gene_ID', 'Seq_ID', 'Type', 'species'])

# simplify header
os.chdir(dir_Seqout)
# os.system("sed -i 's/.*gene=gene:([^\s]+).*/>\\1/' *.fasta")
os.system("sed -i -E 's/.*gene=gene:([[:alnum:]]+).*/>\\1/' *.fasta")

# OrthoFinder
os.chdir(re.sub(r'^(.*\/)[^\/]+\/?$', r'\1', dir_Seqout))
os.system('nohup orthofinder -z -d -a 5 -S blast_nucl -s '+ file_tre + ' -o ' + dir_out + ' -f ' + dir_Seqout + ' &')
#os.system('cp ' + dir_Seqin[:-4] + 'aa/* ' + dir_Seqout[:-4])
