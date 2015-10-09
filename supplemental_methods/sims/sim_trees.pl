#!/n/sw/fasrcsw/apps/Core/perl/5.10.1-fasrc04/bin/perl

use strict;
use warnings;

my $nreps=1000;

#get rates
my $rate = shift;

my $ratedir = sprintf("%.5f", $rate);
system("mkdir -p run2/$ratedir");
for (my $i=0; $i<=$nreps; $i++) {
	#simulate data
	my $simcommand = "java -XX:+UseSerialGC -Xmx8000M -jar ./jprime-0.3.4.jar GuestTreeGen -minper 0 -min 4 -maxper 10000 -max 10000 ultra.final.nwk $rate $rate 0 ./run2/$ratedir/sim.$i";
	system($simcommand);
}

