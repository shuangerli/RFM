# RBP Functional Target Predictor using Multiple Peak Callers (RFM)
Clone the repo: `git clone https://github.com/shuangerli/bio_practicum.git`.

Make the script executable (if necessary): `chmod 755 run.sh`.

Run the script with: `run [argument]`, `./run [argument]` or `slurm run [argument]`.

Please use `run -h` to see instructions.

For an example of the pipeline process, try `run -example`.

Note:

1. If running example, and the downloaded files are broken, please download them manually following these commands:

```
mkdir rep1 rep2
wget -O rep1/reads.R1.fastq.gz https://www.encodeproject.org/files/ENCFF173ODI/@@download/ENCFF173ODI.fastq.gz
wget -O rep1/reads.R2.fastq.gz https://www.encodeproject.org/files/ENCFF507EGU/@@download/ENCFF507EGU.fastq.gz
wget -O rep2/reads.R1.fastq.gz https://www.encodeproject.org/files/ENCFF578XQC/@@download/ENCFF578XQC.fastq.gz
wget -O rep2/reads.R2.fastq.gz https://www.encodeproject.org/files/ENCFF735QHY/@@download/ENCFF735QHY.fastq.gz
wget -O annotation.gtf.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.primary_assembly.annotation.gtf.gz
gunzip annotation.gtf.gz
wget -O ref.GRCh38.fa.gz ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
gunzip ref.GRCh38.fa.gz
wget -O clipper.bed.gz https://www.encodeproject.org/files/ENCFF862EKD/@@download/ENCFF862EKD.bed.gz
wget -O control.bam https://www.encodeproject.org/files/ENCFF608WNY/@@download/ENCFF608WNY.bam
```
