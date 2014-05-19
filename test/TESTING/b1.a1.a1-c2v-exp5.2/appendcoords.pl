#!/usr/bin/perl

open (FILE, "coord.in");
@array1 = grep(/\n/i,<FILE>);
close (FILE);

open (FILE2, ">>attempted.coords");
print FILE2 "\n";
print FILE2 @array1;
close (FILE2);

