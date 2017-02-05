#!/usr/bin/perl

my $filenameOut=$ARGV[2];
my $infile1=$ARGV[0];
my $infile2=$ARGV[1];

my $seqnum=0;
open my $OUTFILE,">$filenameOut";
open my $IN1,"<$infile1";
open my $IN2,"<$infile2";
while(<$IN1>)
{
	if($_=~/^@([\S\s]*)/)
	{
		$_=">".$seqnum.":".$1;
		print $OUTFILE $_;
		$_=<$IN1>;
		print $OUTFILE $_;
		$_=<$IN1>;
		$_=<$IN1>;
		print $OUTFILE "#\n";
		$seqnum++;
	 }
	while(<$IN2>)
	{
		if($_=~/^@([\S\s]*)/)
		{
			$_=">".$seqnum.":".$1;
			print $OUTFILE $_;
			$_=<$IN2>;
			print $OUTFILE $_;
			$_=<$IN2>;
			$_=<$IN2>;
			print $OUTFILE "#\n";
			$seqnum++;
			last;
	 	}
	}#while: <$IN2>
}# while: <$In1>

close $OUTFILE;
close $In1;
close $In2;
