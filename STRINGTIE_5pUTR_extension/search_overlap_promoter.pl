#!/Users/lorenziha/miniconda3/envs/TK_85/bin/perl
use strict;

my $usage = "$0 -i <intergenic transcript gtf> -p <promoter gff> [-w # , default=100 -f <feature>, def = transcript]\n\n";
my %arg = @ARGV;
die $usage unless $arg{-i} && $arg{-p};

my $WINDOW = $arg{-w} || 100;
my $FEATURE = $arg{-f} || 'transcript';

# data structure
# chromosome -> strand -> ([e5, e3], ...)
#

# Open promoter file

my (%cp);
open (my $fhi, "<$arg{-p}") || die "ERROR, I cannot open $arg{-p}: $!\n\n";
while(<$fhi>){
	chomp;
	next if m/^#/;

	my ($chr, undef, $feat, $e5, $e3, $sc, $strand, undef, $comm) = split /\t/;
	my $cp_id = $1 if $comm =~ /ID=(\S+?);/;
	
	push @{$cp{$chr}->{$strand}}, [$e5, $e3, $sc];

}
close $fhi;

# Open transcript gtf# 
my ($t5, $t3); # coord flanking start of transcript
open (my $fhi, "<$arg{-i}") || die "ERROR, I cannot open $arg{-i}: $!\n\n";
while(<$fhi>){
	chomp;
	next if m/^#/;
	

	my ($chr, undef, $feat, $e5, $e3, $sc, $strand, undef, $comm) = split /\t/;
	next unless $feat eq $FEATURE;
	next if $chr eq 'Mito' || $chr eq 'M'; # Exclude mitochondrial chromosome
	my $entry = $_;
	
	# Expand 
	if ($strand eq '.'){

		# try + strand
		($t5, $t3) = ($e5 - $WINDOW, $e5 + $WINDOW);
		# Look for overlapping promoters on same strand
		my @array = @{ $cp{$chr}->{'+'} };
		foreach my $e (@array){
			if ( $e->[1] <= $t5 && $e->[0] >= $t3 ){
				# overlap
				print "$entry @{$e}\n"; 
			}
		}

		# try - strand
		($t5, $t3) = ($e3 - $WINDOW, $e3 + $WINDOW);
		# Look for overlapping promoters on same strand
		my @array = @{ $cp{$chr}->{'-'} };
		foreach my $e (@array){
			if ( $e->[1] <= $t5 && $e->[0] >= $t3 ){
				# overlap
				print "$entry @{$e}\n"; 
			}
		}

	} else {
		if ($strand eq '+'){ 
			($t5, $t3) = ($e5 - $WINDOW, $e5 + $WINDOW);
		} else {
			($t5, $t3) = ($e3 - $WINDOW, $e3 + $WINDOW);
		}	
		
		
		# Look for overlapping promoters on same strand
		my @array = @{ $cp{$chr}->{$strand} };
		foreach my $e (@array){
			my $E5 = $e->[0];
			my $E3 = $e->[1];
			#print "$chr $strand <$e5 - $e3> ( $t5 - $t3 ) >> [ $E5 - $E3 ]\n";

			if ( $E3 >= $t5 && $E5 <= $t3 ){
				# overlap
				print "$entry @{$e}\n"; 
			}
		}
	}
}
close $fhi;

