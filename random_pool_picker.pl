#!/usr/bin/perl
use strict; use warnings;


my @alive = (23..24);
my @dead = (25..26);

assign(@alive);
assign(@dead);



sub assign
{
my @A1;
my @A2;
my $size;
	foreach (@_){
		$size = (@_ / 2);
		my $range = 2;
		my $random = int(rand($range));
		if ($random == 0) {
			if (@A1 < $size){
				push @A1, $_;
			}
			else{
				push @A2, $_;
			}
		}
		if ($random == 1) {
			if (@A2 < $size){
				push @A2, $_;
			}
			else{
				push @A1, $_;
				
			}
		}
	}

print "@A1\n@A2\n"; 
}




