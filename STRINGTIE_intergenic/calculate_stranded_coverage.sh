#!/bin/zsh

# It requires samtools and bedtools
# conda activate TK_85

PREFIX=$1
BAM=$2

# read paired, mapped in proper pair, and first in pair
# NOT read unmapped, NOT mate unmapped, NOT not primary alignment, NOT secondary alignment
# flip strand

echo Running samtools view -@ 8 -f 67 -F 2316 -b $BAM
time samtools view -@ 8 -f 67 -F 2316 -b $BAM > ${PREFIX}_first_filtered.bam

echo Running bedtools genomecov -d -strand + -ibam ${PREFIX}_first_filtered.ba
time bedtools genomecov -d -strand + -ibam ${PREFIX}_first_filtered.bam \
        > ${PREFIX}_minus_3p.cov

cat ${PREFIX}_minus_3p.cov | \
	perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (-1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_minus_3p.bedgraph

echo Running bedtools genomecov -d -strand - -ibam ${PREFIX}_first_filtered.bam
time bedtools genomecov -d -strand - -ibam ${PREFIX}_first_filtered.bam \
        > ${PREFIX}_plus_3p.cov

cat ${PREFIX}_plus_3p.cov | \
        perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_plus_3p.bedgraph


# read paired, mapped in proper pair, and second in pair
# NOT read unmapped, NOT mate unmapped, NOT not primary alignment, NOT secondary alignment

echo Running samtools view -@ 10 -f 131 -F 2316 -b $BAM
time samtools view -@ 10 -f 131 -F 2316 -b $BAM > ${PREFIX}_second_filtered.bam

echo Running bedtools genomecov -d -strand + -ibam ${PREFIX}_second_filtered.bam
time bedtools genomecov -d -strand + -ibam ${PREFIX}_second_filtered.bam \
        > ${PREFIX}_plus_5p.cov

cat ${PREFIX}_plus_5p.cov | \
        perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_plus_5p.bedgraph

echo Runing bedtools genomecov -d -strand - -ibam  ${PREFIX}_second_filtered.bam
time bedtools genomecov -d -strand - -ibam  ${PREFIX}_second_filtered.bam \
        > ${PREFIX}_minus_5p.cov

cat ${PREFIX}_minus_5p.cov | \
        perl -lane '$e5=$F[1]-1; $F[2] = $F[2] * (-1) ;print "$F[0] $e5 $F[1] $F[2]"' > ${PREFIX}_minus_5p.bedgraph


