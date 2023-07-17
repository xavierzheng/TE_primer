#!/bin/bash

set -e -u -o pipefail

less potential_amplicon_summary.txt |awk '$2==1{print $1"_"}' > specific_set.txt

grep -F -f specific_set.txt primer_target_site.bed > specific_primer_target_site.bed

cat specific_set.txt |sed 's/_//g'|awk '{print $1"\t"}' |grep -F -f - potential_amplicon.bed > specific_potential_amplicon.bed
