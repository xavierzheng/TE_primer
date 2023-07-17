#!/bin/bash

# we will add TE ID to primer set name in specific_potential_amplicon.bed and specific_primer_target_site.bed
# modify LTR12_R to LTRXXXXXlot12_R, and XXXXX is TE ID

set -e -u -o pipefail

prefix=$1

cd ${prefix}

mv specific_potential_amplicon.bed specific_potential_amplicon.bed.temp
cat specific_potential_amplicon.bed.temp | awk -v prefix=${prefix} 'BEGIN{FS="\t";OFS="\t"}{gsub("LTR", "lot", $4); print $1, $2, $3, prefix""$4, $5, $6}' > specific_potential_amplicon.bed

mv specific_primer_target_site.bed specific_primer_target_site.bed.temp
cat specific_primer_target_site.bed.temp | awk -v prefix=${prefix} 'BEGIN{FS="\t";OFS="\t"}{gsub("LTR", "lot", $4); print $1, $2, $3, prefix""$4, $5, $6}' > specific_primer_target_site.bed
