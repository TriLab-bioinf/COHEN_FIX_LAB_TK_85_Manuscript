#!/bin/bash

set -e

# Activate conda environment
# conda activate tk_85

# Upload config file with paths to:
# BAM files having the following format: ${SAMPLE_PREFIX}_*.sorted.dedup.bam
# Reference genome file Saccharomyces_cerevisiae.R64-1-1.106.fa
# Annotation file saccharomyces_cerevisiae_R64-4-1_20230830.gff
# Annotation file with annotations of promoter clusters from YeasTSS DB ScerYPDconsensusClusters.gff
# Annotated UTR file Nagalakshmi_2008_5_3_UTRs_V64.bed

source config.txt

# Aux functions
function repeat {
	for idx in $(seq -s " " ${1}); do
		echo -n "="
	done
	echo
}

function my_message {
	len=${#1}
	repeat $len
	echo $1
	repeat $len
	echo
}

# Set path to correct openjdk version (17.0.10 or higher)
export JAVA_HOME=/opt/anaconda3/envs/tk_85/lib/jvm

# Prediction of novel intergenic transcripts

./gff2gtf.pl ${ANNOTATION} > saccharomyces_cerevisiae_R64-4-1_20230830_with_repeats.gtf

# Extract LTR coords from SGD gff file
grep LTR ${ANNOTATION} | perl -ne '@x=split /\t/;$x[3]--;$name=$1 if $x[8]=~m/Name=(\S+?);/ ;print "$x[0]\t$x[3]\t$x[4]\t$name\t.\t$x[6]\n"' > LTR.bed

#for i in {1..12} {31..42}; do
for i in {1..2}; do
	echo Processing sample AY${i}

	# Run stranded stringtie to predict transcripts	
	if [[ ! -e AY${i}.gtf ]]; then
		my_message "Run stranded stringtie to predict transcripts"
		stringtie ${BAM}/AY${i}_*.sorted.dedup.bam --rf -o AY${i}.gtf -G ${ANNOT_GTF} -p 4
	fi

	# Extract transcript IDs that overlap with LTR sequences
	if [[ ! -e AY${i}.LTR.ids ]]; then
		my_message "Extract transcript IDs that overlap with LTR sequences"
		bedtools intersect -a AY${i}.gtf -b LTR.bed -wa|perl -ne 'print "$1\n" if m/gene_id\s"(\S+?)";/' > AY${i}.LTR.ids
	fi

	# Keep only transcripts that do not overlap with LTRs
	if [[ ! -e AY${i}.no_LTR.gtf ]]; then
		my_message "Keep only transcripts that do not overlap with LTRs"
		grep -vwf AY${i}.LTR.ids AY${i}.gtf > AY${i}.no_LTR.gtf
	fi

	# Remove potential novel genes that look like variants of annotated genes
	if [[ ! -e AY${i}.no_LTR.NO_ref_genes.gtf ]]; then
		my_message "Remove potential novel genes that look like variants of annotated genes"
		grep ref_gene_id AY${i}.no_LTR.gtf|perl -ne 'print "$1\n" if m/gene_id\s"(\S+?)";/' > AY${i}.no_LTR.ref_genes.ids
		cat AY${i}.no_LTR.gtf | \
			perl -e 'open I, "AY1.no_LTR.ref_genes.ids"; 
					while(<I>){
						chomp;
						$h{$_}++
					}; 
					close I; 
					while(<>){
						$geneid="none"; 
						$geneid=$1 if m/gene_id\s"(\S+?)";/ ; 
						print "$_" unless $h{$geneid} > 0
					}' > AY${i}.no_LTR.NO_ref_genes.gtf
	fi

	# Generate stranded bigwig file 
	if [[ ! -e AY${i}.wig ]]; then
		my_message "Generate stranded bigwig file"
		igvtools count --strands first -z 7 -w 15  ${BAM}/AY${i}_*.sorted.dedup.bam AY${i}.wig ${REFERENCE}
	fi

	# Substract depth on rev strand from fwd strand
	if [[ ! -e AY${i}.both.wig ]]; then
		my_message "Substract depth on rev strand from fwd strand"
		cat AY${i}.wig | perl -ne 'chomp;@x=split /\t/; if($x[0] =~ m/^\d/ ){$cov = $x[2]-$x[1];print "$x[0]\t$cov\n"}else{print "$x[0]\n"}' > AY${i}.both.wig
	fi

	# Convert wig file to bigwig format
	if [[ ! -e AY${i}.both.bw ]]; then
		my_message "Convert wig file to bigwig format"
		./wigToBigWig AY${i}.both.wig ${CHR_SIZE} AY${i}.both.bw
	fi

	# Remove no_LTR.NO_ref_genes genes that still overlap with annotated genes (taking into account the strand)
	if [[ ! -e AY${i}.no_overlap.no_LTR.NO_ref_genes.gtf ]]; then
		my_message "Remove no_LTR.NO_ref_genes genes that still overlap with annotated genes (taking into account the strand)"
		bedtools intersect -a AY${i}.no_LTR.NO_ref_genes.gtf -b ${ANNOT_GTF} -wa -v -s > AY${i}.no_overlap.no_LTR.NO_ref_genes.gtf
	fi

done

# Merge all overlapping novel transcripts across all samples and save it into aymerge.no_LTR.NO_ref_genes.gtf
my_message "Merge all overlapping novel transcripts across all samples"
stringtie --merge AY*.no_overlap.no_LTR.NO_ref_genes.gtf  > aymerge.no_LTR.NO_ref_genes.gtf

# Identify predicted transcripts with a predicted promoter within a window of 100bp up/downstream from the 5' end
my_message "Identify predicted transcripts with a predicted promoter within a window of 100bp up/downstream from the 5end"
./search_overlap_promoter.pl -i aymerge.no_LTR.NO_ref_genes.gtf -p ${TSS} -w 100 | perl -ne 'print "$1\n" if m/gene_id\s"(MSTRG.\d+?)";/ ' > aymerge.no_LTR.NO_ref_genes_with_PROMOTERS_100bp_window.txt

# Detect predicted merged intergenic transcripts that overlap with predicted UTRs by Nagalakshmi 2008 v64 by at least 10% (-F 0.1) for both -a and -b (-r) and that have the same strand (-s)
my_message "Detect predicted merged intergenic transcripts that overlap with predicted UTRs by Nagalakshmi 2008 v64"
bedtools intersect -a aymerge.no_LTR.NO_ref_genes.gtf -b ${UTR} -wa -s -f 0.1 -r -wb |grep -v exon| perl -ne '@x=split /\t/;print "$1\t$x[12]\n" if m/(MSTRG.\d+)/' > aymerge.no_LTR.NO_ref_genes_Overlap_Nagalakshmi_2008_5_3_UTRs_V64.txt

my_message "Pipeline complete!"  
