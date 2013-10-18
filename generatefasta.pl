#generate new fasta outputs from txt files for input into the Decombinator
use warnings;
use strict;

my $file;
my $path    = 'J:\LAST CD_4';
my $pathtwo = 'J:\LAST CD_4';
my @LSfiles = qw(O2b4-7_F4-3.DNASeq.txt
);
my @filelist;

opendir( DIR, $path ) or die "can't opendir $path: $!";

while ( defined( $file = readdir(DIR) ) ) {
	next if $file =~ /^\.\.?$/;    # skip . and ..
	for my $i ( 0 .. $#LSfiles ) {
		my $tomatch = $LSfiles[$i];
		if ( $file =~ /$tomatch/ ) {
			push( @filelist, $file );    #read in all files from directory
			print("$file\n");
		}
	}
}
closedir(DIR);

for my $t ( 0 .. $#filelist ) {
	my $now     = time;
	my $counter = 0;
	print("Currently on file: $filelist[$t]\n");
	my @data
	  ;   #data is an array of arrayrefs to 1x2 arrays of AA sequence and value;
	my @listofrandomarrays;
	my $currentfile = $filelist[$t];
	my $toreplace   = '.txt';
	$currentfile =~ s/$toreplace//;
	my $currentoutputfiller = "$currentfile" . "_seq" . '.fastq';
	my $currentoutput       = "$pathtwo/$currentoutputfiller";
	open( INPUT,  "$path/$filelist[$t]" );
	open( OUTPUT, ">$currentoutput" );
	my $line = <INPUT>;
	print "$line";

	while ( my $line = <INPUT> ) {
		$line =~ s/\n$//;
		$counter++;
		my @completedata;
		my @currentarray = split( /\t/, $line );

		#			if ($currentarray[$#currentarray] == 0 )
		#			{
		#my $revdna = revdnacomp( $currentarray[4] );
		push( @completedata, $currentarray[3] );

		#print("Currently adding $currentarray[3] to my file\n");
		push( @completedata, $currentarray[4]  );

		#print("Currently adding $currentarray[6] to my file\n");
		push( @completedata, $currentarray[7] );

		#print("Currently adding $currentarray[5] to my file\n");
		my $toprint =
		  to_fastq( $completedata[0], $completedata[1], $completedata[2] );
		print OUTPUT "$toprint";

		#			}
	}
	close(INPUT);
	close(OUTPUT);

	#	open( OUTPUT, ">$currentoutput" );
	#	for my $i (0..$#data)
	#	{
	#		my $toprint = to_fastq($data[$i]->[0],$data[$i]->[1], $data[$i]->[2]);
	#		print OUTPUT "$toprint";
	#	}
	#	close(OUTPUT);
	print("Completed $currentfile with $counter lines");
	$now = time - $now;
	printf(
		"\n\nTotal running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 )
	);

}
my $counter = 0;

sub to_fastq {
	my ( $seqName, $seq, $quality ) = @_;
	$quality =~ s/\s$//;
	$seqName =~ s/\s$//;
	my @qualityarray = split( /,/, $quality );
	my $toprint = "";
	my $seqlength = length($seq);
	for my $i ( 0 .. $seqlength-1) {
		

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

sub revdnacomp {
	my $dna     = shift;
	my $revcomp = reverse($dna);

	$revcomp =~ tr/ACGTacgt/TGCAtgca/;

	return $revcomp;
}
