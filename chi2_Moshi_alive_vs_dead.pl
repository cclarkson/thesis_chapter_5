#!/usr/bin/perl

##########################################################################
##########    											                               #######
########## CHI-SQUARE PROBABILITY CALCULATOR FOR ALL MOSHI (5 vs 3)#######
########## WITHIN vs BETWEEN ALIVE/DEAD COMPARISONS                #######
########## INPUT = All_Moshi_Igor_2L/2R etc_SUBSAMPLED.sync        #######
########## OUTPUT = probabilities based on 2x2 allele count        #######
##########	(biallelic - major/minor alive/dead	                   #######
##########	four comarisons:			                                 ####### 							 
##########	A1vA2,D1vD2,A1vD1,A2vD2		                             #######
##########	pools assigned using random_pool_picker.pl             #######
##########	A1 = 1, 4, 5, 6, 10			                   		         #######
########## 	A2 = 2, 3, 7, 8, 9 			                               #######
##########	D1 = 11, 13, 14				                                 #######
##########	D2 = 12, 15, 16				                                 #######
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

print OUTPUT "contig\tposition\tfisherA1A2\tfisherD1D2\tfisherA1D1\tfisherA2D2\n";


while(my $line = <INPUT>){
	chomp($line);
	my @alleles = split(/\t/, $line);
	my $contig = $alleles[0];
	my $position = $alleles [1];
	shift @alleles for (0..2);
	my @A1 = @alleles[0,3,4,5,9];
	my @A2 = @alleles[1,2,6,7,8];
	my @D1 = @alleles[10,12,13];
	my @D2 = @alleles[11,14,15];

	my $A11; my $A12; my $A13; my $A14;
	my $A21; my $A22; my $A23; my $A24;
	my $D11; my $D12; my $D13; my $D14;
	my $D21; my $D22; my $D23; my $D24;
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
		foreach(@D1){
			my @D1alleles = split(/:/, $_);
			$D11 += $D1alleles[0];
			$D12 += $D1alleles[1];
			$D13 += $D1alleles[2];
			$D14 += $D1alleles[3];
		}
		foreach(@D2){
			my @D2alleles = split(/:/, $_);
			$D21 += $D2alleles[0];
			$D22 += $D2alleles[1];
			$D23 += $D2alleles[2];
			$D24 += $D2alleles[3];
		}
	
		my $total1 = $A11 + $A21 + $D11 + $D21;
		my $total2 = $A12 + $A22 + $D12 + $D22;
		my $total3 = $A13 + $A23 + $D13 + $D23;
		my $total4 = $A14 + $A24 + $D14 + $D24;

		push (my @totalcount, $total1, $total2, $total3, $total4);
		my @sortedtotal = sort {$a <=> $b } @totalcount;
		if ($sortedtotal[1] == 0){ #only use biallelics

			push( my @A1_all_count, $A11 , $A12 , $A13, $A14);
                        push( my @A2_all_count, $A21 , $A22 , $A23, $A24);
                        push( my @D1_all_count, $D11 , $D12 , $D13, $D14);
                        push( my @D2_all_count, $D21 , $D22 , $D23, $D24);

			my @A1_non_zero;
			my @A2_non_zero;
			my @D1_non_zero;
			my @D2_non_zero;
			

			if ($A1_all_count[0] != 0 || $A2_all_count[0] != 0 || $D1_all_count[0] != 0 || $D2_all_count[0] != 0) {
		        chomp($A1_all_count[0]);
        		chomp($A2_all_count[0]);
        		chomp($D1_all_count[0]);
        		chomp($D2_all_count[0]);
        		push (@A1_non_zero, $A1_all_count[0]);
        		push (@A2_non_zero, $A2_all_count[0]);
			push (@D1_non_zero, $D1_all_count[0]);
			push (@D2_non_zero, $D2_all_count[0]);
			}
			if ($A1_all_count[1] != 0 || $A2_all_count[1] != 0 || $D1_all_count[1] != 0 || $D2_all_count[1] != 0) {
                        chomp($A1_all_count[1]);
                        chomp($A2_all_count[1]);
                        chomp($D1_all_count[1]);                                                                                      
                        chomp($D2_all_count[1]);                                                                                      
                        push (@A1_non_zero, $A1_all_count[1]);                                                                        
                        push (@A2_non_zero, $A2_all_count[1]);
                        push (@D1_non_zero, $D1_all_count[1]);
                        push (@D2_non_zero, $D2_all_count[1]);                                                                        
                        }			
			if ($A1_all_count[2] != 0 || $A2_all_count[2] != 0 || $D1_all_count[2] != 0 || $D2_all_count[2] != 0) {
                        chomp($A1_all_count[2]);
                        chomp($A2_all_count[2]);
                        chomp($D1_all_count[2]);                                                                                      
                        chomp($D2_all_count[2]);                                                                                      
                        push (@A1_non_zero, $A1_all_count[2]);                                                                        
                        push (@A2_non_zero, $A2_all_count[2]);
                        push (@D1_non_zero, $D1_all_count[2]);
                        push (@D2_non_zero, $D2_all_count[2]);                                                                       
                        }
			if ($A1_all_count[3] != 0 || $A2_all_count[3] != 0 || $D1_all_count[3] != 0 || $D2_all_count[3] != 0) {
			chomp($A1_all_count[3]);
			chomp($A2_all_count[3]);
			chomp($D1_all_count[3]);
			chomp($D2_all_count[3]);
			push (@A1_non_zero, $A1_all_count[3]);
			push (@A2_non_zero, $A2_all_count[3]);
			push (@D1_non_zero, $D1_all_count[3]);
			push (@D2_non_zero, $D2_all_count[3]);
			}

			if(@A1_non_zero < 2){
				push(@A1_non_zero, 0);
			} 
			if(@A2_non_zero < 2){
				push(@A2_non_zero, 0);
			}
			if(@D1_non_zero < 2){
				push(@D1_non_zero, 0);
			}
			if(@D2_non_zero < 2){
				push(@D2_non_zero, 0);
			}



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
			

                        my @fixDD1;
                        my @fixDD2;
                        if ($D1_non_zero[0] <= $D1_non_zero[1]){
                        @fixDD1 = @D1_non_zero;
                        @fixDD2 = @D2_non_zero; 
                        }
                        if ($D1_non_zero[0] > $D1_non_zero[1]){
			@fixDD1 = sort {$a <=> $b} @D1_non_zero; 
			@fixDD2 = ($D2_non_zero[1],$D2_non_zero[0])
			}
			my $inp1DD = $fixDD1[0];
			my $inp2DD = $fixDD1[1];
			my $inp3DD = $fixDD2[0];
			my $inp4DD = $fixDD2[1];
			my @finalDD = ([$inp1DD, $inp2DD], [$inp3DD, $inp4DD]);
			my $DDchi = new Statistics::ChisqIndep;
			$DDchi->load_data(\@finalDD);
			my $fisherD1D2 = $DDchi->{p_value};
			


                        my @fix1AD1;
                        my @fix1AD2;
                       if ($A1_non_zero[0] <= $A1_non_zero[1]){
                        @fix1AD1 = @A1_non_zero;
                        @fix1AD2 = @D1_non_zero;
                        }
                        if ($A1_non_zero[0] > $A1_non_zero[1]){
                        @fix1AD1 = sort {$a <=> $b} @A1_non_zero;
                        @fix1AD2 = ($D1_non_zero[1],$D1_non_zero[0]);
                        }
			my $inp1AD1 = $fix1AD1[0];
			my $inp2AD1 = $fix1AD1[1];
			my $inp3AD1 = $fix1AD2[0];
			my $inp4AD1 = $fix1AD2[1];
			my @finalAD1 = ([$inp1AD1, $inp2AD1], [$inp3AD1, $inp4AD1]);
			my $AD1chi = new Statistics::ChisqIndep;
			$AD1chi->load_data(\@finalAD1);
			my $fisherA1D1 = $AD1chi->{p_value};
		
			my @fix2AD1;
                        my @fix2AD2;
                        if ($A2_non_zero[0] <= $A2_non_zero[1]){
                        @fix2AD1 = @A2_non_zero;
                        @fix2AD2 = @D2_non_zero;
                        }
                        if ($A2_non_zero[0] > $A2_non_zero[1]){
                        @fix2AD1 = sort {$a <=> $b} @A2_non_zero;
                        @fix2AD2 = ($D2_non_zero[1],$D2_non_zero[0]);
                        }
			my $inp1AD2 = $fix2AD1[0];
			my $inp2AD2 = $fix2AD1[1];
			my $inp3AD2 = $fix2AD2[0];
			my $inp4AD2 = $fix2AD2[1];
			my @finalAD2 = ([$inp1AD2, $inp2AD2], [$inp3AD2, $inp4AD2]);
			my $AD2chi = new Statistics::ChisqIndep;
			$AD2chi->load_data(\@finalAD2);
			my $fisherA2D2 = $AD2chi->{p_value};
			
				print OUTPUT "$contig\t$position\t$fisherA1A2\t$fisherD1D2\t$fisherA1D1\t$fisherA2D2\n";
				

	}
}





