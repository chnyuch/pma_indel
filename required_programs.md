### Download software
alias bcftools='/vol/storage/bcftools'

#### Git

```bash
# for git clone
conda install anaconda::git
```

#### Blast

```bash
conda install bioconda::blast
```

#### NCBI-genome-download

```bash
conda install -c bioconda ncbi-genome-download
```

#### NCBI dataset

https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

```bash
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
```

#### OrthoFinder

```bash
conda install bioconda::orthofinder
```

#### AGAT

https://github.com/NBISweden/AGAT

```bash
conda create -n agat -c bioconda agat
```

#### Prank

```bash
conda install bioconda::prank
```

#### indelMap

https://github.com/acg-team/indelMaP

```bash
cd /vol/storage/software
git clone https://github.com/acg-team/indelMaP.git
python /vol/storage/software/indelMaP/src/indelMaP_ASR.py
```

#### seqkit

https://bioinf.shenwei.me/seqkit/

```bash
conda install bioconda::seqkit
```

#### pal2nal

https://www.bork.embl.de/pal2nal/

```bash
conda config --env --add channels anaconda
conda install anaconda::gcc_linux-64
conda install bioconda::pal2nal
```

#### bedtools

```bash
conda install bioconda::bedtools
```

#### paml

https://github.com/abacus-gene/paml

http://abacus.gene.ucl.ac.uk/software/#paml-for-unixlinux

```bash
wget https://github.com/abacus-gene/paml/releases/download/4.10.7/paml-4.10.7-linux-X86_64.tgz
tar zxvf paml-4.10.7-linux-X86_64.tgz
```

#### 

### 
