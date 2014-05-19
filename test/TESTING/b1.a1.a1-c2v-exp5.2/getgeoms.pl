#!/usr/bin/perl
#
use strict;
use warnings;

# This script grabs the coordinates for a geometry from the
#  surfgen input

# getgeoms.pl -g [geometry] -a [atoms]

my $geom = $ARGV[1];
my $atoms= $ARGV[3];

open(FILE1,"refgeom") or die "No file.";
my @file = grep(/\n/i, <FILE1>);
close (FILE1);

my $i=0;
my $j=0;
my @g1;
for ($i = ($geom - 1)*($atoms+1); $i <= ($geom -1)*($atoms+1) + $atoms; $i++){
        $g1[$j] = $file[$i];
	$j = $j+1;
}
print @g1;
