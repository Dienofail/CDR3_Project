#script to classify the V_CDR3_J regions of Dr. Yongqing's GBset.xlsx file.
#The output of this script is meant to be fed directly to decombinator_edited_again.py

use warnings;
use strict;
use List::Compare;
use Test::More tests => 1;
my $inputfile = "GBset.txt";
my $bad_file = "n_list.txt";
my @sequence;
my @bad_list;
open( INPUT, "$inputfile" );
my $header = <INPUT>;
while ( my $line = <INPUT> ) {
	$line =~ s/\n$//;
	my @current_line = split( /\t/, $line );
	if ( $#current_line > 1 ) {
		push( @sequence, $current_line[5] );
	}
}
close(INPUT);

my $outputfile   = "GBset_decombinator_input" . ".fastq";
my $name_counter = 0;
my $bad_counter  = 0;
open( OUTPUT, ">$outputfile" );
for my $i ( 0 .. $#sequence ) {
	if ( $sequence[$i] =~ m/\d/ ) {
		next;
	}
	my $name_counter++;
	my $length  = length( $sequence[$i] );
	my @chars   = split( "", $sequence[$i] );
	my $counter = 0;

	for my $i ( 0 .. $#chars ) {
		if ( $chars[$i] =~ m/[ATGC\s]/ ) {
			$counter++;
		}
	}

	if ( $counter != $length ) {
		$bad_counter++;
		push(@bad_list, $sequence[$i]);
		print "$sequence[$i]\n";
		next;
	}
	my $quality = '40';
	for my $i ( 2 .. $length ) {
		$quality .= ',' . "40";
	}
	my @quality_counter = split( /,/, $quality );
	my $seq_name_input = '@' . "$name_counter";

	#is($#quality_counter+1,length($sequence[$i]));
	#print("$bad_counter\n");
	my $fastq_to_print = to_fastq( $seq_name_input, $sequence[$i], $quality );
	print OUTPUT "$fastq_to_print";
}
print("$bad_counter\n");
close(OUTPUT);
my @bad_list_ref;
open (INPUT, "$bad_file");
my $orig_num = 0;
while (my $line = <INPUT>)
{
	$orig_num++;
	$line =~ s/\n$//;
	push (@bad_list_ref, $line);
}
close(INPUT);
my $lc = List::Compare->new(\@bad_list_ref, \@bad_list);

my @intersection = $lc->get_intersection;
my @L_only = $lc->get_Lonly;
my @R_only = $lc->get_Ronly;

print("Orig_Num from decombinator: $orig_num Intersection: $#intersection Ref only: $#L_only New: $#R_only\n");
for my $i (0..$#R_only)
{
	print "$R_only[$i]\n";
}
sub to_fastq {
	my ( $seqName, $seq, $quality ) = @_;
	$quality =~ s/\s$//;
	$seqName =~ s/\s$//;
	my @qualityarray = split( /,/, $quality );
	my $toprint      = "";
	my $seqlength    = length($seq);
	my $counter      = 0;
	for my $i ( 0 .. $seqlength - 1 ) {

		#print("Currently processing $qualityarray[$i]\n");
		my $topush = $qualityarray[$i] + 33;

		#print("Added 32 to my score, generating $topush\n");
		my $topushchar = chr($topush);

		#print("Converted topush to char, generating $topushchar\n");
		$toprint = $toprint . $topushchar;

		#print("I am extending toprint with $toprint\n");
		$counter++;
	}

	# default to 80 characters of sequence per line
	#        my $len = 78;
	#        for my $j (0..77)
	#        {
	#        	$toprint = $toprint . 'A';
	#        }
	my $formatted_seq = "$seqName" . "$counter\n";
	$formatted_seq .= "$seq\n";
	$seqName =~ s/\@/\+/;
	$formatted_seq .= "$seqName" . "$counter\n";
	$formatted_seq .= "$toprint\n";
	return $formatted_seq;
}
