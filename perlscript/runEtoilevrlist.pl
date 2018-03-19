#!/usr/bin/perl -w

########################
# Run ETOILE to select spectra from MILES to measure RV 

#create list of all .ascii to run 
`ls rv*.ascii > rvlist`;
#inicializations 
my $i = 0;
my $name = "";

open (INFILE, 'rvlist');
while(<INFILE>){
chomp;
# to run etoile for all stars (1 step)
#`$name= 'echo rvTer9-0266_802.ascii | sed s/.ascii//' `
$name=substr($_,0,13);
`printf "%d\n%s\n%s\n%d\n" 0 $_ $name 5  | ./etoile`;

print "$_..........OK!\n";
$i=$i+1;
}
close INFILE;
