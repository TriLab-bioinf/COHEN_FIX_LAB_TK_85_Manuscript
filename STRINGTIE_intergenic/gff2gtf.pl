#!/usr/bin/perl
use strict;

# I   SGD chromosome  1   230218  .   .   .   ID=chrI;dbxref=NCBI:BK006935.2;Name=chrI
# I   SGD telomere    1   801 .   -   .   ID=TEL01L;Name=TEL01L;Note=Telomeric region on the left arm of Chromosome I  composed of an X element core sequence  X element combinatorial repeats  and a short terminal stretch of telomeric repeats;display=Telomeric region on the left arm of Chromosome I;dbxref=SGD:S000028862;curie=SGD:S000028862
# I   SGD X_element   337 801 .   -   .   Parent=TEL01L;Name=TEL01L_X_element
# I   SGD X_element_combinatorial_repeat  63  336 .   -   .   Parent=TEL01L;Name=TEL01L_X_element_combinatorial_repeat
# I   SGD telomeric_repeat    1   62  .   -   .   Parent=TEL01L;Name=TEL01L_telomeric_repeat_1
# I   SGD telomere    1   801 .   -   .   ID=TEL01L_telomere;Name=TEL01L_telomere;Parent=TEL01L
# I   SGD gene    335 649 .   +   .   ID=YAL069W;Name=YAL069W;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious open reading frame  unlikely to encode a functional protein  based on available experimental and comparative sequence data;display=Dubious open reading frame;dbxref=SGD:S000002143;orf_classification=Dubious;curie=SGD:S000002143
# I   SGD CDS 335 649 .   +   0   Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classification=Dubious;protein_id=UniProtKB:O13588
# I   SGD mRNA    335 649 .   +   .   ID=YAL069W_mRNA;Name=YAL069W_mRNA;Parent=YAL069W
# I   SGD gene    538 792 .   +   .   ID=YAL068W-A;Name=YAL068W-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Dubious open reading frame  unlikely to encode a functional protein  based on available experimental and comparative sequence data  identified by gene-trapping  microarray-based expression analysis  and genome-wide homology searching;display=Dubious open reading frame;dbxref=SGD:S000028594;orf_classification=Dubious;curie=SGD:S000028594
# I   SGD CDS 538 792 .   +   0   Parent=YAL068W-A_mRNA;Name=YAL068W-A_CDS;orf_classification=Dubious;protein_id=UniProtKB:Q8TGK7
# I   SGD mRNA    538 792 .   +   .   ID=YAL068W-A_mRNA;Name=YAL068W-A_mRNA;Parent=YAL068W-A
# I   SGD ARS 707 776 .   .   .   ID=ARS102;Name=ARS102;Alias=ARSI-1;Note=Autonomously Replicating Sequence;display=Autonomously Replicating Sequence;dbxref=SGD:S000121252;curie=SGD:S000121252
# I   SGD gene    1807    2169    .   -   .   ID=YAL068C;Name=YAL068C;gene=PAU8;Alias=PAU8,seripauperin PAU8;Ontology_term=GO:0003674,GO:0005575,GO:0030437,SO:0000704;Note=Protein of unknown function  member of the seripauperin multigene family encoded mainly in subtelomeric regions;display=Protein of unknown function;dbxref=SGD:S000002142;orf_classification=Verified;curie=SGD:S000002142
# I   SGD CDS 1807    2169    .   -   0   Parent=YAL068C_mRNA;Name=YAL068C_CDS;orf_classification=Verified;protein_id=UniProtKB:P0CE92
# I   SGD mRNA    1807    2169    .   -   .   ID=YAL068C_mRNA;Name=YAL068C_mRNA;Parent=YAL068C;transcript_id=RefSeq:NM_001181127.1
# I   SGD gene    2480    2707    .   +   .   ID=YAL067W-A;Name=YAL067W-A;Ontology_term=GO:0003674,GO:0005575,GO:0008150,SO:0000704;Note=Putative protein of unknown function  identified by gene-trapping  microarray-based expression analysis  and genome-wide homology searching;display=Putative protein of unknown function;dbxref=SGD:S000028593;orf_classification=Uncharacterized;curie=SGD:S000028593


# Reads gff from stdin
# Outputs gtf to stdout


#   31 ncRNA_gene
#   12 pseudogene
#   24 rRNA_gene
#    6 snRNA_gene
#   77 snoRNA_gene
#  299 tRNA_gene
#    1 telomerase_RNA_gene
#   91 transposable_element_gene

# outputs
# I	SGD	CDS	335	649	.	+	0	transcript_id "YAL069W_mRNA"; gene_id "YAL069W"; gene_name "YAL069W";
# I	SGD	CDS	538	792	.	+	0	transcript_id "YAL068W-A_mRNA"; gene_id "YAL068W-A"; gene_name "YAL068W-A";


my ($chr,$db,$feat,$e5,$e3,$score,$strand,$phase,$com, %types);

my @types = ("ncRNA_gene", "pseudogene", "rRNA_gene", "snRNA_gene", "snoRNA_gene", "tRNA_gene", 
			"telomerase_RNA_gene", "transposable_element_gene", "LTR_retrotransposon", "long_terminal_repeat");

foreach my $type (@types){
	$types{$type}++
} 

while(<>){
	
	# gff header
	if(m/^#/){
		print "$_";
		next;
	}
	
	chomp;
	($chr,$db,$feat,$e5,$e3,$score,$strand,$phase,$com) = split /\t/;
	

	next unless ($types{$feat});

	my ($t_id, $g_id, $g_name) = &get_info_from_com($com);
	my $com_parsed = "gene_id \"$g_id\"; transcript_id \"$t_id\"; gene_name \"$g_name\"";

	print "$chr\t$db\ttranscript\t$e5\t$e3\t$score\t$strand\t$phase\t$com_parsed\n";
	
}

sub get_info_from_com {
	my $com = shift;
	my ($transcript_id, $gene_id, $gene_name) = ("","","");
	my $transcript_id = $1."_mRNA" if $com =~ m/ID=(\S+?);/;
	my $gene_id = $1 if $com =~ m/Name=(\S+?);/;
	my $gene_name = $2 if $com =~ m/(gene|Note)=([\S\s]+?);/;

	return($transcript_id,$gene_id,$gene_name);

}
