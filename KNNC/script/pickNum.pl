#!/usr/bin/perl -w

my $infile=$ARGV[0];
my $copynum=$ARGV[1];
my $seqnum=$ARGV[2];

open IN,"<$infile";
open OUT1,">$copynum";
open OUT2,">$seqnum";

while(<IN>)
{
	if($_=~/Frequency:([\d]*)=>VertexNum:([\d]*)/)
	 {
		 print OUT1 $1."\n";
		 print OUT2 $2."\n";
	 }
}

close IN;
close OUT1;
close OUT2;
