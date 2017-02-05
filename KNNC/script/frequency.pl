#!/usr/bin/perl -w

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my %count;

open IN,"<$infile";
open OUT,">$outfile";

print "$infile\n";
print "$outfile\n";

while(<IN>)
{
   if($_=~/frequency=([\d]*)/)
	{
	   if($1)  
	   {
	        ++$count{$1};
	   } 
	}
}

for $key (sort {$a <=> $b} keys %count)
{
	my $str="Frequency:".$key."=>"."VertexNum:".$count{$key}."\n";
	print OUT $str;
}

print "Finished!\n";
close IN;
close OUT;
