#!/usr/bin/perl -w
use strict;
use feature "switch";

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my %commonFreq;
my %commonDegree;
my $No=0;

open IN,"<$infile";
open OUT,">$outfile";

print "$infile\n";
print "$outfile\n";


my $str;
my $id;
my $freq;
my $degree;
my $maxkmer=0;
$str=<IN>;
if($str=~/k-num:([\d]*)/)
{
	$maxkmer=$1;
}
$str=<IN>;
my $realkmer=0;
if($str=~/kmer:([\d]*)/)
{
	 $realkmer=$1;
}


while(<IN>)
{
   if($_=~/>Id=([\d]*);frequency=([\d]*)/)
	{
	  $id=$1;	
	  $freq=$2;
	   if($freq>20)  
	   {
			my $seq=<IN>;
			my $k=0;
			$seq=<IN>;
			$id=&getId($seq);
			while(<IN>)
			{
	   			if($_=~/Neighbors:([\d]*)/){
					
	          		$commonFreq{$id}=$freq;
					$commonDegree{$id}=$1;
					last;
				}#if
			}#while
	   }#if 
	}#if
}

$str="No.\tId\tFrequency*1000000\tDegree\n";
print OUT $str;
my $frate=0;
my $vrate=0;
$No=1;
for my $key (sort {$a <=> $b} keys %commonFreq)
{
	#my $str=$No."\t".$key."\t".$commonFreq{$key}."\t".$commonDegree{$key}."\n";
	my $frate=($commonFreq{$key}/$maxkmer)*1000000;
	#my $vrate=$count{$key}/$size;
	#printf OUT "%0.5f\t%0.5f\n",$frate,$vrate;
	my $str=$No."\t".$key."\t".$frate."\t".$commonDegree{$key}."\n";
	print OUT $str;
	$No++;
}

close IN;
close OUT;

sub getId
{
	my $seq=shift;
	my @str=split "",$seq;
	my $k=shift;
	my $id=0;
	my $i=0;
	my $k=$#str;
	
	print "k:$k\n";

	foreach my $ch(@str)
	{
		given ($ch)
		{
			when($ch=~/A/) {$id+=0*&power(4,($k-$i-1));}
			when($ch=~/C/) {$id+=1*&power(4,($k-$i-1));}
			when($ch=~/G/) {$id+=2*&power(4,($k-$i-1));}
			when($ch=~/T/) {$id+=3*&power(4,($k-$i-1));}
		}
		$i++;
	}
	return $id;
}

sub power
{
	my $a=shift;
	my $b=shift;
	my $pow=1;

	for(my $i=0;$i<$b;$i++)
	{
		$pow=$pow * $a;
	}
	return $pow;
}
