#!/bin/bash
# ReadSeeker Evaluation Download Data


# Download Reference Sequences
mkdir -p references/raw
mkdir -p references/human
mkdir -p references/sars2
mkdir -p references/gammaherpes
mkdir -p references/mtuberkolosis
mkdir -p references/ecoli
mkdir -p references/mouse


echo "Downloading References"

# human GRCh38
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gbff.gz -P references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_000001405.40_GRCh38.p14_genomic.gbff.gz references/human/human

# mouse genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gbff.gz -p references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_000001635.27_GRCm39_genomic.gbff.gz references/mouse/mouse

# M. tuberculosis NC_000962.3
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gbff.gz -P references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_000195955.2_ASM19595v2_genomic.gbff.gz references/mtuberkolosis/mtuberkolosis


# E. coli NC_000913.3
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz -P references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_000005845.2_ASM584v2_genomic.gbff.gz references/ecoli/ecoli

# Epstein Barr Virus GCF_002402265
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/402/265/GCF_002402265.1_ASM240226v1/GCF_002402265.1_ASM240226v1_genomic.gbff.gz -P references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_002402265.1_ASM240226v1_genomic.gbff.gz references/gammaherpes/gammaherpes

# Sars-CoV-2 NC_045512.2
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gbff.gz -P references/raw
python scripts/pythonCDS_BLAST.py makecds references/raw/GCF_009858895.2_ASM985889v3_genomic.gbff.gz references/sars2/sars2


echo "Generating Blastdb"

# Make Tblastn db
makeblastdb -in references/mtuberkolosis/mtuberkolosis_genome.fasta -dbtype nucl -out references/mtuberkolosis/mtuberkolosis_nucl_blastdb
makeblastdb -in references/gammaherpes/gammaherpes_genome.fasta -dbtype nucl -out references/gammaherpes/gammaherpes_nucl_blastdb
makeblastdb -in references/ecoli/ecoli_genome.fasta -dbtype nucl -out references/ecoli/ecoli_nucl_blastdb



echo "load uniprot sequences"
wget https://zenodo.org/records/12699064/files/uniprot_seq.zip?download=1 -O references/uniprot.zip
unzip -d references/ references/uniprot.zip

echo "Blasting Protein Sequences"

# Extract uniprot sequences
zcat references/uniprot_seq/uniprotkb_taxonomy_id_1762_2023_11_23.fasta.gz > references/mtuberkolosis/mtuberkolosis_uniprotkb_taxonomy_id_1762_2023_11_23.fasta
zcat references/uniprot_seq/uniprotkb_taxonomy_id_548681_2023_12_15.fasta.gz > references/gammaherpes/gammaherpes_uniprotkb_taxonomy_id_548681_2023_12_15.fasta
zcat references/uniprot_seq/uniprotkb_taxonomy_id_561_2024_01_12.fasta.gz > references/ecoli/ecoli_uniprotkb_taxonomy_id_561_2024_01_12.fasta

# Blast Protein Sequences to Reference Genomes
tblastn -query references/mtuberkolosis/mtuberkolosis_uniprotkb_taxonomy_id_1762_2023_11_23.fasta -db references/mtuberkolosis/mtuberkolosis_nucl_blastdb -outfmt 6 -num_threads 55 -out references/mtuberkolosis/mtuberkolosis_blastoutput.tsv
tblastn -query references/gammaherpes/gammaherpes_uniprotkb_taxonomy_id_548681_2023_12_15.fasta -db references/gammaherpes/gammaherpes_nucl_blastdb -outfmt 6 -num_threads 55 -out references/gammaherpes/gammaherpes_blastoutput.tsv
tblastn -query references/ecoli/ecoli_uniprotkb_taxonomy_id_561_2024_01_12.fasta -db references/ecoli/ecoli_nucl_blastdb -outfmt 6 -num_threads 55 -out references/ecoli/ecoli_blastoutput.tsv

# Generate CDS files from Blast result

echo "Generating CDS files"
python scripts/pythonCDS_BLAST.py blast2bed references/mtuberkolosis/mtuberkolosis_cds.bed  references/mtuberkolosis/mtuberkolosis_blastoutput.tsv references/mtuberkolosis/mtuberkolosis_uniprotkb_taxonomy_id_1762_2023_11_23.fasta references/mtuberkolosis/mtuberkolosis_cds_uniprotkb_taxonomy_id_1762_2023_11_23_joined.bed
python scripts/pythonCDS_BLAST.py blast2bed references/gammaherpes/gammaherpes_cds.bed references/gammaherpes/gammaherpes_blastoutput.tsv references/gammaherpes/gammaherpes_uniprotkb_taxonomy_id_548681_2023_12_15.fasta references/gammaherpes/gammaherpes_cds_uniprotkb_taxonomy_id_548681_2023_12_15_joined.bed
python scripts/pythonCDS_BLAST.py blast2bed references/ecoli/ecoli_cds.bed references/ecoli/ecoli_blastoutput.tsv references/ecoli/ecoli_uniprotkb_taxonomy_id_561_2024_01_12.fasta references/ecoli/ecoli_cds_uniprotkb_taxonomy_id_561_2024_01_12_joined.bed








# Download Samples
mkdir samples
echo downloading samples

# Human 1
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/082/ERR10492982/ERR10492982_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/082/ERR10492982/ERR10492982_2.fastq.gz -P samples/

# Human 2
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/041/ERR10493241/ERR10493241_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR104/041/ERR10493241/ERR10493241_2.fastq.gz -P samples/

#Human 3
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/072/ERR10509672/ERR10509672_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/072/ERR10509672/ERR10509672_2.fastq.gz -P samples/

# M. tuberculosis 1
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/022/SRR21820122/SRR21820122_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/022/SRR21820122/SRR21820122_2.fastq.gz -P samples/

# M. tuberculosis 2
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/024/SRR21820124/SRR21820124_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/024/SRR21820124/SRR21820124_2.fastq.gz -P samples/

# M. tuberculosis 3
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/055/SRR21864655/SRR21864655_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR218/055/SRR21864655/SRR21864655_2.fastq.gz -P samples/

# E. coli
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR226/087/SRR22674487/SRR22674487_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR226/087/SRR22674487/SRR22674487_2.fastq.gz -P samples/

# Epstein Barr Virus
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/008/ERR2024408/ERR2024408_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR202/008/ERR2024408/ERR2024408_2.fastq.gz -P samples/

# Sars-CoV-2 1
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/059/ERR10913059/ERR10913059_1.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/059/ERR10913059/ERR10913059_2.fastq.gz -P samples/

# Sars-CoV-2 2
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/061/ERR10913061/ERR10913061_2.fastq.gz -P samples/
wget -nc http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/061/ERR10913061/ERR10913061_1.fastq.gz -P samples/

# Mouse
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR317/DRR317657/DRR317657_1.fastq.gz -P samples/
wget http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR317/DRR317657/DRR317657_2.fastq.gz -P samples/

# Trim mouse sample to 2 million reads (1mio Read 1, 1mio Read 2)
zcat DRR317657_2.fastq.gz | head -n 20000000 | gzip > DRR317657_sub5m_2.fastq.gz
zcat DRR317657_1.fastq.gz | head -n 20000000 | gzip > DRR317657_sub5m_1.fastq.gz



zcat samples/ERR3021870_1.fastq.gz | head -n 4000000 | gzip > ERR3021870sub1m_1.fastq.gz 
zcat samples/ERR3021870_2.fastq.gz | head -n 4000000 | gzip > ERR3021870sub1m_2.fastq.gz



# build singularity image
#singularity build singularity/dnabert.sif singularity/dnabert.def 
