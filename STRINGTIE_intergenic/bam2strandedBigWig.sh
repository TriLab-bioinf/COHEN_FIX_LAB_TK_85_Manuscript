#!/bin/bash
#

PREFIX=$1
WD=/Users/lorenziha/Documents/DKBIOCORE_LOCAL/TK_85/STRINGTIE
BAM_PATH=/Users/lorenziha/Documents/DKBIOCORE_LOCAL/TK_85/data/BAM
GENOME=Saccharomyces_cerevisiae.R64-1-1.106.fa

# Generate stranded big wig file

BAM=$(ls ${BAM_PATH}/${PREFIX}_S*_L*.sorted.dedup.bam)

echo BAM file = ${BAM}

time igvtools count --strands first -z 7 -w 25  ${BAM} ${PREFIX}.wig ${GENOME}

cat ${PREFIX}.wig | perl -ne 'chomp;@x=split /\t/; if($x[0] =~ m/^\d/ ){$cov = $x[2]-$x[1];print "$x[0]\t$cov\n"}else{print "$x[0]\n"}' > ${PREFIX}.both.wig

${WD}/wigToBigWig ${PREFIX}.both.wig ${WD}/chromosome_sizes.txt ${PREFIX}.both.bw



