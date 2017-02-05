#!/usr/bin/perl

my $filenum=$#ARGV;
my $filenameOut=$ARGV[$filenum];
my $seqnum=0;
my $num=0;
open my $OUTFILE,">$filenameOut";
for($num=0;$num<$filenum;$num++)
{
	my $filename=shift ;
	open my $INFILE,"<$filename";
	while(<$INFILE>)
	{
		if($_=~/^@([\S\s]*)/)
		{
		   $_=">".$seqnum.":".$1;
		   print $OUTFILE $_;
		   $_=<$INFILE>;
		   print $OUTFILE $_;
		   $_=<$INFILE>;
		   $_=<$INFILE>;
		   print $OUTFILE "#\n";
		   $seqnum++;
		   next;
		}
		if($_=~/^>([\S\s]*)/)
		{
			$_=">".$seqnum."|".$1;
			print $OUTFILE $_;
			while(<$INFILE>)
			{
				if($_=~/^>([\S\s]*)/)
				  {
				     print $OUTFILE "#\n";
					 $seqnum++;
					 $_=">".$seqnum."|".$1;
					 print $OUTFILE $_;
				  }
			    else
				  {
					print $OUTFILE $_;
				  }
			 }
				print $OUTFILE "#";
				$seqnum++;
		}
	}# while
	close $INFILE;
}# for 

close $OUTFILE;
