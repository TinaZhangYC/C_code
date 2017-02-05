#!/usr/bin/perl

$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

my $seqNum=0;

open $FILEA, "< $filenameA";
open $FILEB, "< $filenameB";

open $OUTFILE, "> $filenameOut";

while(<$FILEA>) {
	$_=">".$seqNum."\n";
	print $OUTFILE $_;
	$seqNum++;
	$_ = <$FILEA>;
	print $OUTFILE $_;
	$_ = <$FILEA>;
	$_ = <$FILEA>;
	print $OUTFILE "#\n";
}
while(<$FILEB>){
	$_=">".$seqNum."\n";
	print $OUTFILE $_;
	$seqNum++;
	$_ = <$FILEB>;
	print $OUTFILE $_; 
	$_ = <$FILEB>;
	$_ = <$FILEB>;
	print $OUTFILE "#\n";
}
close $FILEA;
close $FILEB;
vlose $OUTFILE;
