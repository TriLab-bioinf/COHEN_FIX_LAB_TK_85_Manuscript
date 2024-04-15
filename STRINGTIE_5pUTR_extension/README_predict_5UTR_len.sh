#!/bin/bash

# Activate conda environment
# conda activate tk_85

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

#### Predict 5'UTR lengths for each protein-coding gene for each individual strand 


BAM=/Users/lorenziha/Documents/LOCAL_PROJECTS/ORNA_COHEN_FIX/TK_85/data/BAM

my_message "Make stranded coverage in bed and bedgraph formats"
for i in {1..12} {31..42}; do 
	./calculate_stranded_coverage_v2.sh AY${i} ${BAM}/AY${i}_*.sorted.dedup.bam
done

exit

my_message "Concatenate bed files for plus and minus strands"
cat AY1_plus.bed > all.plus.bed
cat AY1_minus.bed > all.minus.bed

for i in {2..12} {31..42}; do
	cut -f 3 AY${i}_plus.bed | paste all.plus.bed - > tmp1
	mv tmp1 all.plus.bed

	cut -f 3 AY${i}_minus.bed | paste all.minus.bed - > tmp2
	mv tmp2 all.minus.bed
done

my_message "Predict 5pUTR lengths using the following cutoffs: sequencing gap = 0 and minimim sequencing depth = 1"
./look_for_5p_extension.pl -g ./Saccharomyces_cerevisiae.R64-1-1.106.minus.gtf -d all.minus.bed -w 200 -n 24 -D 1 -G 0 > all.minus.200bp.D1.G0.5p_extension
./look_for_5p_extension.pl -g ./Saccharomyces_cerevisiae.R64-1-1.106.plus.gtf  -d all.plus.bed  -w 200 -n 24 -D 1 -G 0 > all.plus.200bp.D1.G0.5p_extension

my_message "format 5p_extension files as bed files for display in IGV"
./5p_extension_to_bed.pl -i all.plus.200bp.D1.G0.5p_extension -s plus -w 200 > all.plus.200bp.D1.G0.5p_extension.bed
./5p_extension_to_bed.pl -i all.minus.200bp.D1.G0.5p_extension -s minus -w 200 > all.minus.200bp.D1.G0.5p_extension.bed

my_message "Filter WT annotations from 5p_extension.bed files"
grep 'AY9\|AY10\|AY11\|AY12\|AY39\|AY40\|AY41\|AY42' all.plus.200bp.D1.G0.5p_extension.bed > all.plus.200bp.D1.G0.WT.5p_extension.bed
grep 'AY9\|AY10\|AY11\|AY12\|AY39\|AY40\|AY41\|AY42' all.minus.200bp.D1.G0.5p_extension.bed > all.minus.200bp.D1.G0.WT.5p_extension.bed

