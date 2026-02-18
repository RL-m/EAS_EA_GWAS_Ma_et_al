source ~/.bashrc
source-conda
source /packages/miniconda3/etc/profile.d/conda.sh
conda activate ldsc
inFile=$1
outFile=$2
munge_sumstats="/packages/ldsc/munge_sumstats.py"
snpList="/LDSC/reference/w_eas_hm3.snplist"
${munge_sumstats} \
--sumstats ${inFile} \
--signed-sumstat BETA,0 --snp SNP --p P --N-col N \
--merge-alleles ${snpList} --chunksize 500000 \
--out ${outFile}