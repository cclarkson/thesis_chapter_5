#!/usr/bin/perl
##########################################################################
##########    											                               #######
########## CHI-SQUARE PROBABILITY CALCULATOR FOR MA v Mer (6 v 1)  #######
########## INPUT = All_Moshi_Igor_2L/2R etc_SUBSAMPLED.sync        #######
########## OUTPUT = probabilities based on 2x2 allele count        #######
##########	(biallelic - two comparisons:			                     #######							 
##########	A1vA2		                                               #######
##########	MA = 11,12,13,14,15,16		                             #######
########## 	Merus = 23 			                                       #######
##########			by Chris Clarkson	                                 #######
##########			csc@liv.ac.uk		                                   #######
##########################################################################





use strict; 
use warnings;
use Statistics::ChisqIndep;
use POSIX;


print "Pleez inputz the name ov output file\n\n";

chomp (my $title = <STDIN>);


my $input  = $ARGV[0];
open (INPUT, $input)||die "cannot open input";
open (OUTPUT, ">$title".".txt")||die "aaaaaaaaaaaaarrrrrrrrrrgh!";
my $counter = 0;

print OUTPUT "contig\tposition\tchi2MoshiDeadvMerus\n";


while(my $line = <INPUT>){
	chomp($line);
	my @alleles = split(/\t/, $line);
	my $contig = $alleles[0];
	my $position = $alleles [1];
	shift @alleles for (0..2);
	my @A1 = @alleles[10,11,12,13,14,15];
	my @A2 = @alleles[22];

	my $A11; my $A12; my $A13; my $A14;
	my $A21; my $A22; my $A23; my $A24;

		foreach(@A1){
                	my @A1alleles = split(/:/, $_);
                	$A11 += $A1alleles[0];
                	$A12 += $A1alleles[1];
                	$A13 += $A1alleles[2];
                	$A14 += $A1alleles[3];
        	}
		foreach(@A2){
                	my @A2alleles = split(/:/, $_);
                	$A21 += $A2alleles[0];
                	$A22 += $A2alleles[1];
                	$A23 += $A2alleles[2];
                	$A24 += $A2alleles[3];
        	}
	
		my $total1 = $A11 + $A21;
		my $total2 = $A12 + $A22;

		#print "$total1\t$total2\t$total3\t$total4\n";  #get true ATCG counts across all pools to get rid of >biallelics
		push (my @totalcount, $total1, $total2);
		my @sortedtotal = sort {$a <=> $b } @totalcount;
		if ($sortedtotal[1] == 0){ #only use biallelics

			push( my @A1_all_count, $A11 , $A12 , $A13, $A14);
                        push( my @A2_all_count, $A21 , $A22 , $A23, $A24);

			my @A1_non_zero;
			my @A2_non_zero;
			

			if ($A1_all_count[0] != 0 || $A2_all_count[0] != 0 ) {
		        chomp($A1_all_count[0]);
        		chomp($A2_all_count[0]);
        		push (@A1_non_zero, $A1_all_count[0]);
        		push (@A2_non_zero, $A2_all_count[0]);
			}
			if ($A1_all_count[1] != 0 || $A2_all_count[1] != 0 ) {
                        chomp($A1_all_count[1]);
                        chomp($A2_all_count[1]);                                                                                      
                        push (@A1_non_zero, $A1_all_count[1]);                                                                        
                        push (@A2_non_zero, $A2_all_count[1]);                                                                        
                        }			
			if ($A1_all_count[2] != 0 || $A2_all_count[2] != 0 ) {
                        chomp($A1_all_count[2]);
                        chomp($A2_all_count[2]);                                                                                      
                        push (@A1_non_zero, $A1_all_count[2]);                                                                        
                        push (@A2_non_zero, $A2_all_count[2]); 
                        }
			if ($A1_all_count[3] != 0 || $A2_all_count[3] != 0 ) {
			chomp($A1_all_count[3]);
			chomp($A2_all_count[3]);
			push (@A1_non_zero, $A1_all_count[3]);
			push (@A2_non_zero, $A2_all_count[3]);
			}

			if(@A1_non_zero < 2){
				push(@A1_non_zero, 0);
			} 
			if(@A2_non_zero < 2){
				push(@A2_non_zero, 0);
			}
			#print "@A1_non_zero\t@A2_non_zero\t@D1_non_zero\t@D2_non_zero\n";


			my @fixAA1;
			my @fixAA2;
			if ($A1_non_zero[0] <= $A1_non_zero[1]){
		        @fixAA1 = @A1_non_zero;
		        @fixAA2 = @A2_non_zero;
			}
			if ($A1_non_zero[0] > $A1_non_zero[1]){
		        @fixAA1 = sort {$a <=> $b} @A1_non_zero;
		        @fixAA2 = ($A2_non_zero[1],$A2_non_zero[0]);
			}
			my $inp1AA = $fixAA1[0];
			my $inp2AA = $fixAA1[1];
			my $inp3AA = $fixAA2[0];
			my $inp4AA = $fixAA2[1];
			my @finalAA = ([$inp1AA, $inp2AA], [$inp3AA, $inp4AA]);
			my $AAchi = new Statistics::ChisqIndep;
			$AAchi->load_data(\@finalAA);
			my $fisherA1A2 = $AAchi->{p_value};

				print OUTPUT "$contig\t$position\t$fisherA1A2\n";

	}
}





