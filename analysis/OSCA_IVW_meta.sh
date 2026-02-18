OSCA="/packages/osca-0.46.1-linux-x86_64/osca-0.46.1"
fileList=$1
outFile=$2
${OSCA} \
--gwas-flist ${fileList} \
--meta \
--out ${outFile}