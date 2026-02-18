# calculate disease score
h5adFile=$1
gsFile=$2
covFile=$3
mainAdd=$4
scdrs="/packages/scDRS/bin/scdrs"
${scdrs} compute-score \
    --h5ad-file ${h5adFile} \
    --h5ad-species human \
    --gs-file ${gsFile} \
    --gs-species human \
    --cov-file ${covFile} \
    --flag-filter-data True \
    --flag-raw-count True \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True \
    --out-folder ${mainAdd}/

# calculate cell-type association
h5adFile=$1
scoreFile=$2 #input full_score
outAdd=$3
group=$4
scdrs="/packages/scDRS/bin/scdrs"
${scdrs} perform-downstream \
        --h5ad-file ${h5adFile} \
        --score-file ${scoreFile} \
        --out-folder ${outAdd}/ \
        --group-analysis ${group} \
        --flag-filter-data True \
        --flag-raw-count True