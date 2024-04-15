#!/Users/lorenziha/miniconda3/envs/TK_85/bin/perl
use strict;

## Example input gtf from stringtie:

#I   StringTie   transcript  7235    9016    1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; reference_id "YAL067C_mRNA"; ref_gene_id "YAL067C"; ref_gene_name "SEO1"; cov "1.394876"; FPKM "0.528678"; TPM "0.577279";
#I   StringTie   exon    7235    9016    1000    -   .   gene_id "STRG.1"; transcript_id "STRG.1.1"; exon_number "1"; reference_id "YAL067C_mRNA"; ref_gene_id "YAL067C"; ref_gene_name "SEO1"; cov "1.394876";
#I   StringTie   transcript  11565   11951   1000    -   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; reference_id "YAL065C_mRNA"; ref_gene_id "YAL065C"; cov "23.608534"; FPKM "8.947971"; TPM "9.770553";
#I   StringTie   exon    11565   11951   1000    -   .   gene_id "STRG.2"; transcript_id "STRG.2.1"; exon_number "1"; reference_id "YAL065C_mRNA"; ref_gene_id "YAL065C"; cov "23.608534";
#I   StringTie   transcript  12046   12426   1000    +   .   gene_id "STRG.3"; transcript_id "STRG.3.1"; reference_id "YAL064W-B_mRNA"; ref_gene_id "YAL064W-B"; cov "9.800912"; FPKM "3.714686"; TPM "4.056174";
#I   StringTie   exon    12046   12426   1000    +   .   gene_id "STRG.3"; transcript_id "STRG.3.1"; exon_number "1"; reference_id "YAL064W-B_mRNA"; ref_gene_id "YAL064W-B"; cov "9.800912";
#I   StringTie   transcript  13363   13743   1000    -   .   gene_id "STRG.4"; transcript_id "STRG.4.1"; reference_id "YAL064C-A_mRNA"; ref_gene_id "YAL064C-A"; ref_gene_name "TDA8"; cov "15.632972"; FPKM "5.925119"; TPM "6.469812";
#I   StringTie   exon    13363   13743   1000    -   .   gene_id "STRG.4"; transcript_id "STRG.4.1"; exon_number "1"; reference_id "YAL064C-A_mRNA"; ref_gene_id "YAL064C-A"; ref_gene_name "TDA8"; cov "15.632972";

my $usage = "$0 -g <gtf file>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-g};

my (%t);
open (my $fhi,"<$arg{-g}") || die "ERROR, I cannot open $arg{-g}: $!\n\n";
while(<$fhi>){
    chomp;
    my ($chr, $db, $feat, $e5, $e3, $sc, $strand, $phase, $comm) = split /\t/;
    next unless $feat eq "transcript";

    my ($gene_id, $ref_gene_id, $gene_name) = &parse_comm($comm);
    if ($e5 < $t{$gene_id}->{e5} || !exists($t{$gene_id}->{e5}) ){$t{$gene_id}->{e5} = $e5;};
    if ($e3 > $t{$gene_id}->{e3} || !exists($t{$gene_id}->{e3})){$t{$gene_id}->{e3} = $e3;};
    $t{$gene_id}->{strand} = $strand;
    $t{$gene_id}->{name} = $gene_name || "none";
    $t{$gene_id}->{chr} = $chr;
    $t{$gene_id}->{ref_gene_id} = $ref_gene_id || "none";
}
close $fhi;

# Print output
foreach my $gene_id (keys %t){
    my $e5 = $t{$gene_id}->{e5};
    my $e3 = $t{$gene_id}->{e3};
    my $strand = $t{$gene_id}->{strand};
    my $ref_gene_name = $t{$gene_id}->{name};
    my $chr = $t{$gene_id}->{chr};
    my $ref_gene_id = $t{$gene_id}->{ref_gene_id};
    if ($ref_gene_id ne 'none'){
        print "$chr\tStringTie_merged\ttranscript\t$e5\t$e3\t.\t$strand\t.\tgene_id \"$gene_id\"; ref_gene_id \"$ref_gene_id\"; ref_gene_name \"$ref_gene_name\";\n";
    }
}

sub parse_comm {
    my $comm = shift;
    $comm =~ s/"//g;
    my %annot = split(m/;?\s/,$comm);

    return($annot{gene_id}, $annot{ref_gene_id}, $annot{ref_gene_name});
}
