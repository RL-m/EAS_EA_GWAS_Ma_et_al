# SNP	A1	A2	N	beta	SE
# Notice the order of sFile: POP1=EUR, POP2=EAS #
popcorn=/packages/Popcorn/popcorn
score=/packages/POPCORN/EUR_EAS_all_gen_eff.cscore
#score=/packages/POPCORN/EUR_EAS_all_gen_imp.cscore
sFile1=$1
sFile2=$2
outFile=$3
popcorn fit -v 1 --cfile ${score} --sfile1 ${sFile1} --sfile2 ${sFile2} --maf 0.01 --gen_effect ${outFile}