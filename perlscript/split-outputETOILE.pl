#!/usr/bin/perl -w

use Cwd;

my $homedir = `echo \$HOME`;
my $workdir = getcwd();

my $i = 1;
my $j = 1;
my $k = 1;

open (INFILE,"<atmosphericParameters.dat");
my @data = <INFILE>;
my $lines=@data;
close (INFILE);

for($i=1;$i<$lines;$i=$i+1){
  if($i==($j-1)*152 + 1){
    print "$data[$i-1]";
    chomp($data[$i-1]);
    open (BUFF,">$data[$i-1]");
    for($k=$k=152*($j-1) + 1;$k<152*($j-1) + 151;$k=$k+1){
      print BUFF "$data[$k]";
    }
    close(BUFF);
    $j=$j+1;
  }
}

