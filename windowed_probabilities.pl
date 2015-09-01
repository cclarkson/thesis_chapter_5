#!/usr/bin/perl

#################################################################### 
####	creates windowed geomean/median for Arabiensis GWAS	      ####
####	INPUT = All_moshi_Igor_N_SUBSAMPLED_probabilities       	####
####							                                             	####
####	Chris Clarkson csc@liv.ac.uk				                      ####
####################################################################

use strict; 
use warnings;
use Statistics::Descriptive;

print "\n\nhereth enter thee name ov output\n\n";

chomp(my $title = <STDIN>);

print "\n\nto proceed thou must now define window size\n\n";

chomp(my $winsize = <STDIN>);


my $input = $ARGV[0];

open (INPUT, $input) || die "cor, me input has gone all wonkz";
open (OUTPUT, ">$title".".txt") || die "halp! output down! medic!";

my $counter = 0;
my $start;
my $end;
my $contig;
my @winAA;
my @winDD;
my @winAD1;
my @winAD2;

print OUTPUT "contig\tstart\tend\tmid\tAAmean\tAAvar\tAAmin\tAAmed\tAAgeo\tDDmean\tDDvar\tDDmin\tDDmed\tDDgeo\tAD1mean\tAD1var\tAD1min\tAD1med\tAD1geo\tAD2mean\tAD2var\tAD2min\tAD2med\tAD2geo\n";


while(my $line = <INPUT>){
	if ($line !~ /^c/){
		$counter ++;
		if ($counter < $winsize ) {
			my @split = split(/\t/ ,$line);
			if ($counter == 1) {
				$start = $split[1];
			}
			chomp($split[5]);
			push (@winAA, $split[2]);		 
			push (@winDD, $split[3]);	 
			push (@winAD1, $split[4]);		 
			push (@winAD2, $split[5]);
			
		}
		if ($counter == $winsize ) {
			my @split = split(/\t/ , $line);
			$contig = $split[0];
			$end = $split[1];
			chomp($split[5]);
			push (@winAA, $split[2]);
			push (@winDD, $split[3]);
			push (@winAD1, $split[4]);
			push (@winAD2, $split[5]);
			(my $AAmean, my $AAvar, my $AAmin, my $AAmed, my $AAgeo) = &mean(@winAA);
			(my $DDmean, my $DDvar, my $DDmin, my $DDmed, my $DDgeo) = &mean(@winDD);
			(my $AD1mean, my $AD1var, my $AD1min, my $AD1med, my $AD1geo) = &mean(@winAD1);
			(my $AD2mean, my $AD2var, my $AD2min, my $AD2med, my $AD2geo) = &mean(@winAD2);
			my $mid = ($end + $start)/2;
			print OUTPUT "$contig\t$start\t$end\t$mid\t$AAmean\t$AAvar\t$AAmin\t$AAmed\t$AAgeo\t$DDmean\t$DDvar\t$DDmin\t$DDmed\t$DDgeo\t$AD1mean\t$AD1var\t$AD1min\t$AD1med\t$AD1geo\t$AD2mean\t$AD2var\t$AD2min\t$AD2med\t$AD2geo\n";
			undef @winAA; undef @winDD; undef @winAD1; undef @winAD2;
			$counter = 0;
		}
	}
}


sub mean
{
my $stat = Statistics::Descriptive::Sparse->new();
$stat ->add_data(@_);
my $submean = $stat->mean();
my $subvar = $stat->variance();
my $submin = $stat->min();
my $fullstat = Statistics::Descriptive::Full->new();
$fullstat ->add_data(@_);
my $submedian = $fullstat->median();
my $subgeo = $fullstat->geometric_mean();
return ($submean, $subvar, $submin, $submedian, $subgeo);
} 



close OUTPUT;
close INPUT;

