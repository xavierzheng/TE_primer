#!/bin/bash

set -e -u -o pipefail

module purge
module load blast/2.13.0

INPUT=${1}
OUTPUT=${2}

blastn -query ${INPUT} -db /nfs/project1/index/BoA_ctm_fix_V2.2.fa \
        -task blastn-short \
        -evalue 100 \
        -outfmt "7 qseqid sseqid qlen slen sstrand pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
	-out ${OUTPUT}

