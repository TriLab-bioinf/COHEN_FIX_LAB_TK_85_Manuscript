#!/bin/zsh

PREFIX=$1
BAM=$2

echo '=============================================================='
echo Processing $BAM
echo '=============================================================='

echo Calculating coverage: bedtools genomecov -ibam $BAM -d -du -strand + -split \> ${PREFIX}_minus.bed

time bedtools genomecov -ibam $BAM -d -du -strand + -split > ${PREFIX}_minus.bed
cat ${PREFIX}_minus.bed | perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (-1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_minus.bedgraph

echo Calculating coverage: bedtools genomecov -ibam $BAM -d -du -strand - -split \> ${PREFIX}_plus.bed

time bedtools genomecov -ibam $BAM -d -du -strand - -split > ${PREFIX}_plus.bed
cat ${PREFIX}_plus.bed | perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_plus.bedgraph

echo
echo '============================== DONE =============================='
echo
