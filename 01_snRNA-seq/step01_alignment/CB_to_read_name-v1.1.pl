#!/usr/bin/perl
use strict;
use warnings;

my $usage = <<_EOUSAGE_;
#########################################################################################
#
# Usage:
#       perl $0 fastq1 out.fastq.gz
#
# Note:
#      add cell barcode in read name(fastq) to read (fastq)
#
# Description:  
#
# History:
#  20250131 v1
#  20250131 v1.1 match the quality length to the barcode base length.
#
#########################################################################################
_EOUSAGE_
    ;

##################
#  at least 2 paramiter
die $usage unless (scalar @ARGV == 2);

#for my $i (@ARGV)
# {
# unless (-s $i)
#     {
#       die "Error: the file $i does not exist!\n";
#
#     }
#}

if ($ARGV[0] =~ /.gz$/) {
open(IN1, "gunzip -c $ARGV[0] |") || die "can not open pipe to $ARGV[0]";
}
else {
open(IN1, $ARGV[0]) || die "can’t open $ARGV[0]";
}


open OUT, "| gzip >$ARGV[1]" or die $!;

while ( my $id = <IN1> ) {

    # read 4 lines in seq1!
    my $seq = <IN1>;
	chomp($seq);
	my $plus = <IN1>;
	my $qual = <IN1>;
	chomp($qual);
    
	#get cell barcode and umi
	$id =~ /\S+_([ATCGN]+)_([ATCGN]+)\s/;
	my $bc = $1.$2;
	my $size = length($bc);
	my $bc_qual = 'I' x $size;
	#add to reads
	$seq = $bc.$seq."\n";
	$qual = $bc_qual.$qual."\n";
	
	#output
    print OUT $id;
	print OUT $seq;
	print OUT $plus;
	print OUT $qual;
}
close IN1;
close OUT;
