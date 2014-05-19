#!/usr/bin/perl

# Check if HD_OLD exists
if (!-e "./HD_OLD"){
	die "No directory.";
}

# Check if iteration file exists, if not create it.
if (-e "./hd_iter"){
	my @i;
	open( FILE, "./hd_iter");
	@i = grep(/\n/i,<FILE>);
	close (FILE);
	$k = @i[0];
} else {
	my $k=1;
        open( FILE, ">>./hd_iter");
	print FILE "1\n";
	close (FILE);
}

system("cp hd.data2 ./HD_OLD/hd.data$k");
system("rm ./hd_iter");

my $j = $k + 1;
open( FILE, ">>./hd_iter");
print FILE "$j\n";
close(FILE);

print "Finished.\n";
