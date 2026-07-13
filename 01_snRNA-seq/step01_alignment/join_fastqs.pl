#!/usr/bin/perl
use strict;
use warnings;

my $usage = <<_EOUSAGE_;
#########################################################################################
#
# Usage:
#       perl $0 fastq1 fastq2 out.fastq.gz
#
# Note:
#     merge sequence for 2 fastq gz input, such as add cell barcode(fastq) to read (fastq)
#     output to stout
#
# Description:  
#
# History:
#
#########################################################################################
_EOUSAGE_
    ;

##################
#  at least 2 paramiter
die $usage unless (scalar @ARGV == 3);

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

if ($ARGV[1] =~ /.gz$/) {
open(IN2, "gunzip -c $ARGV[1] |") || die "can not open pipe to $ARGV[1]";
}
else {
open(IN2, $ARGV[1]) || die "can not open $ARGV[1]";
}

open OUT, "| gzip >$ARGV[2]" or die $!;

while ( my $id1 = <IN1> ) {

    # read 4 lines in seq1!
    my $seq1 = <IN1>;
	chomp($seq1);
	my $plus1 = <IN1>;
	my $qual_1 = <IN1>;
	chomp($qual_1);

    # read 4 lines in seq2!
	my $id2 = <IN2>;
    my $seq2 = <IN2>;
	my $plus2 = <IN2>;
	my $qual_2 = <IN2>;
    # merge
	my $seq = $seq1.$seq2;
	my $qual = $qual_1.$qual_2;
	
	#output
    print OUT $id1;
	print OUT $seq;
	print OUT $plus1;
	print OUT $qual;
}
close IN1;
close IN2;
close OUT;
