#!/usr/bin/perl -w

#run ETOILE list 

#create list of all .ascii to run 
`ls *.ascii > list`;
#inicializations 
open(NET, 'net-only.tab');
@net=<NET>;
close(NET);

my $i = 0;
my $name = "";
my $nameplot="";

open (INFILE, 'list');
while(<INFILE>){
chomp;
# to run etoile for all stars (1 step)
$net[$i]=$net[$i]*1;
$name=substr($_,0,13);
`printf "%d\n%s\n%s\n%d\n%d\n%d\n" 0 $_ $name 5  $net[$i] 0 | ./etoile`;

$nameplot= "spec" . $name . ".tab" ;
`mv spec2plot.tab $nameplot`;

print "$_..........OK!\n";
$i=$i+1;
}
close INFILE;
