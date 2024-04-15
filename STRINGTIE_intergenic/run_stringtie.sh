#!/bin/bash

PREFIX=$1
WD=/Users/lorenziha/Documents/DKBIOCORE_LOCAL/TK_85/STRINGTIE
BAM_PATH=/Users/lorenziha/Documents/DKBIOCORE_LOCAL/TK_85/data/BAM
GENOME=${WD}/Saccharomyces_cerevisiae.R64-1-1.106.fa

# Get BAM file name
BAM=$(ls ${BAM_PATH}/${PREFIX}_S*_L*.sorted.dedup.bam)
echo BAM file = ${BAM}

# Run stranded stringtie
time stringtie ${BAM} --rf -o ${PREFIX}.gtf -G Saccharomyces_cerevisiae.R64-1-1.106.gtf -p 4

# Extract transcript IDs that overlap with LTR sequences
time bedtools intersect -a ${PREFIX}.gtf -b ${WD}/LTR.bed -wa|perl -ne 'print "$1\n" if m/gene_id\s"(\S+?)";/' > ${PREFIX}.LTR.ids

# Keep only transcripts that do not overlap with LTRs
time grep -vwf ${PREFIX}.LTR.ids ${PREFIX}.gtf > ${PREFIX}.no_LTR.gtf

# Remove potential novel genes that look like variants of annotated genes
grep 'ref_gene_id' ${PREFIX}.no_LTR.gtf|perl -ne 'print "$1\n" if m/gene_id\s"(\S+?)";/' > ${PREFIX}.no_LTR.ref_genes.ids
cat ${PREFIX}.no_LTR.gtf | perl -e 'open I, "'${PREFIX}'.no_LTR.ref_genes.ids"; while(<I>){chomp;$h{$_}++}; close I; while(<>){$geneid="none"; $geneid=$1 if m/gene_id\s"(\S+?)";/ ; print "$_" unless $h{$geneid} > 0}' > ${PREFIX}.no_LTR.NO_ref_genes.gtf

# Remove no_LTR.NO_ref_genes genes that still overlap with annotated genes (taking into account the strand)
for p in {1..12} {31..42}; do
bedtools intersect -a AY${p}.no_LTR.NO_ref_genes.gtf -b Saccharomyces_cerevisiae.R64-1-1.106.gtf -wa -v -s > AY${p}.no_overlap.no_LTR.NO_ref_genes.gtf
done

# Rename novel transcripts by sample
❯ for p in {1..12} {31..42}; do
cat AY${p}.no_overlap.no_LTR.NO_ref_genes.gtf |sed 's/STRG/STRG_AY'${p}'/g' >> test.gtf
done

# Merge all overlapping novel transcripts across all samples and save it into aymerge.no_LTR.NO_ref_genes.gtf
stringtie --merge AY*.no_overlap.no_LTR.NO_ref_genes.gtf  > aymerge.no_LTR.NO_ref_genes.gtf

# Generate list of predicted mRNA genes with their respective TPMs (all based on stringtie predictions)
for i in {1..12} {31..42}; do
grep 'mRNA.\+ref_gene_id' AY${i}.no_LTR.gtf |grep "\ttranscript\t" |perl -ne '$g=$1 if m/gene_id\s"(\S+?)";/; $tpm=$1 if m/TPM\s"(\d+.\d*?)";/; print $g."_AY'${i}'\t$tpm\n"' >> mRNA_genes_TPMs.txt
done

grep 'mRNA.\+ref_gene_id' AY${i}.no_LTR.gtf |grep "\ttranscript\t" |perl -ne '$g=$1 if m/gene_id\s"(\S+?)";/; $tpm=$1 if m/TPM\s"(\d+.\d*?)";/; $ref=$1 if m/ref_gene_name\s"(\S+?)";/; print $g."_AY'${i}'\t$tpm\t$ref\n"'|sort -k3 -u

for i in {1..12} {31..42}; do
grep 'mRNA.\+ref_gene_id' AY${i}.no_LTR.gtf |grep "\ttranscript\t" |perl -ne '$g=$1 if m/gene_id\s"(\S+?)";/; $tpm=$1 if m/TPM\s"(\d+.\d*?)";/; $ref=$1 if m/ref_gene_name\s"(\S+?)";/; print "$ref\t$tpm\n"'|sort -k1 -u > AY${i}.mRNA.tpm
done


# Concatenate tpm files based on reference gene name
cat AY1.mRNA.tpm | xlookup.pl -d AY2.mRNA.tpm -dri 1 | xlookup.pl -d AY3.mRNA.tpm -dri 1 | xlookup.pl -d AY4.mRNA.tpm -dri 1 | xlookup.pl -d AY5.mRNA.tpm -dri 1 | xlookup.pl -d AY6.mRNA.tpm -dri 1 | xlookup.pl -d AY7.mRNA.tpm -dri 1 | xlookup.pl -d AY8.mRNA.tpm -dri 1 | xlookup.pl -d AY9.mRNA.tpm -dri 1 | xlookup.pl -d AY10.mRNA.tpm -dri 1 | xlookup.pl -d AY11.mRNA.tpm -dri 1 | xlookup.pl -d AY12.mRNA.tpm -dri 1 | xlookup.pl -d AY31.mRNA.tpm -dri 1 | xlookup.pl -d AY32.mRNA.tpm -dri 1 | xlookup.pl -d AY33.mRNA.tpm -dri 1 | xlookup.pl -d AY34.mRNA.tpm -dri 1 | xlookup.pl -d AY35.mRNA.tpm -dri 1 | xlookup.pl -d AY36.mRNA.tpm -dri 1 | xlookup.pl -d AY37.mRNA.tpm -dri 1 | xlookup.pl -d AY38.mRNA.tpm -dri 1 | xlookup.pl -d AY39.mRNA.tpm -dri 1 | xlookup.pl -d AY40.mRNA.tpm -dri 1 | xlookup.pl -d AY41.mRNA.tpm -dri 1 | xlookup.pl -d AY42.mRNA.tpm -dri 1 > ALL.mRNA.tpm


