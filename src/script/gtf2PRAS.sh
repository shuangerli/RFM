## original programmer: Jianan Lin (jianan.lin@jax.org)
## Last edit: Jan 28th, 2019
## Contact by email for issues.
## download gtfToGenePred Utility from http://hgdownload.soe.ucsc.edu/admin/exe/

## $1 is the gtf file downloaded from UCSC or other resources; $2 is the transcript ID and gene name file, especially for those gtf files which don't include gene name information.
## gtf to genePred format
gtfToGenePred $1 -allErrors "$1".tmp
## segment genepred files
python extract_genomic_regions.py "$1".tmp
## genePred format to PRAS format
python genepred2PRAS.py $2 "$1".tmp "$1".transcript.PRASanno
python genepred2PRAS.py $2 "$1".tmp.utr5.genepred "$1".utr5.PRASanno
python genepred2PRAS.py $2 "$1".tmp.cds.genepred "$1".cds.PRASanno
python genepred2PRAS.py $2 "$1".tmp.utr3.genepred "$1".utr3.PRASanno

rm "$1".tmp
rm "$1".tmp.utr5.genepred
rm "$1".tmp.cds.genepred
rm "$1".tmp.utr3.genepred
