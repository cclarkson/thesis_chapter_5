#!/usr/bin/perl

use strict; use warnings;


my $input = $ARGV[0]; #input file - sync
my $output = $ARGV[1]; #name of output
my $maf = $ARGV[2]; #maf cut off - e.g. for 5% type 5

open (INPUT, $input) || die "must stay away from the light";
open (OUTPUT, ">$output") || die "it's too late for me, save yourself";

while (my $line = <INPUT>){
	chomp($line);
	my @alleles = split(/\t/, $line);
	shift @alleles for (0..2);
	my $A = 0; my $T = 0; my $G = 0; my $C = 0;
		foreach(@alleles){
			my @freqs = split(/:/, $_);
                        $A += $freqs[0];
                        $T += $freqs[1];
                        $G += $freqs[2];
                        $C += $freqs[3];
                }
	my $total = $A + $T + $G + $C;
	push (my @allfreq, $A, $T, $G, $C);
	my @sortedfreq = sort {$a <=> $b } @allfreq;
		if ($sortedfreq[1] == 0){ #no triallelics
			my $maftest = $total / 100 * $maf;
				if ($sortedfreq[2] >= $maftest){
					print OUTPUT "$line\n";
				}
		}

	
}		
 
close INPUT;
close OUTPUT;

