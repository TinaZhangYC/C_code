#!/usr/bin/perl -w

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my %degree;

open IN,"<$infile";
open OUT,">$outfile";

print "$infile\n";
print "$outfile\n";

while(<IN>)
{
   if($_=~/Neighbors:([\d]*)/)
	{
	   #if($1)  
	   #{
	        ++$degree{$1};
	   #} 
	}
}

for $key (sort {$a<=>$b} keys %degree)
{
	my $str="NeighborDegree".$key."=>"."VertexNum:".$degree{$key}."\n";
	print OUT $str;
}

print "Finished!\n";

close IN;
close OUT;
