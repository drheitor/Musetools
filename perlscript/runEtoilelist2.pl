#!/usr/bin/perl -w

#run ETOILE list 

#create list of all .ascii to run 
`ls *.ascii > list`;
#inicializations 
my $i = 0;
my $name = "";

open (INFILE, 'list');
while(<INFILE>){
chomp;
# to run etoile for all stars (1 step)

$name=substr($_,0,11);
`printf "%d\n%s\n%s\n%d\n%f\n" 0 $_ $name 5 2.5 | ./etoile`;

print "$_..........OK!\n";
$i=$i+1;
}
close INFILE;
