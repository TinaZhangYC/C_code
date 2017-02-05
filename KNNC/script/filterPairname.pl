#!/usr/bin/perl -w

&frequency($ARGV[0],$ARGV[1],$ARGV[2]);


sub frequency{
	my $infile1=shift;#sequence
	my $infile2=shift; #Name
	my $outfile=shift;
	my $outfile1=$outfile."_1.flt";
	my $outfile2=$outfile."_2.flt";
	my $outfile3=$outfile.".sgl";
	my $Id;
	my %reads;
	my %names;

	open IN1,"<$infile1";
	open IN2,"<$infile2";
	print "In file:$infile1\t$infile2\n";
	while(<IN1>)
	{
   		if($_=~/^>([\d]*)/)
		{
				$Id=$1;
				$_=<IN1>;
				chomp;
	        	$reads{$Id}=$_;
		}
	}#while
	while(<IN2>)
	{
		chomp;
		if($_=~/^>([\d]*)/)
		{
			$Id=$1;
			$names{$Id}=$_;
		}
	}
	close IN1;
	close IN2;

	open OUT1,">$outfile1";
	open OUT2,">$outfile2";
	open OUT3,">$outfile3";
	my $id1;
	my $id2;
	my $str1;
	my $str2;
	foreach my $id (sort {$a <=> $b} keys %reads)
	{
		if($id %2 ==0)
		{
			$id1=$id;
			$id2=$id1+1;
			if(exists $reads{$id2} ) #pairs
			{
				$str1=$names{$id1}."\n".$reads{$id1}."\n";
				$str2=$names{$id2}."\n".$reads{$id2}."\n";
				print OUT1 $str1;
				print OUT2 $str2;
			}
			else #single
			{
				$str1=$names{$id1}."\n".$reads{$id1}."\n";
				print OUT3 $str1;
			}
		}
		else
		{
			$id2=$id;
			$id1=$id2-1;
			if(exists $reads{$id1} ) #pairs,but is printed!
			{
				next;
			}
			else #single
			{
				$str2=$names{$id2}."\n".$reads{$id2}."\n";
				print OUT3 $str2;
			}
		}
	}#foreach

	close OUT1;
	close OUT2;
	close OUT3;
	return ;
}#sub
