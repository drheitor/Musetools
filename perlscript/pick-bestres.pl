#!/usr/bin/perl -w

use Cwd;
#use Env qw(CODE);
#use Switch;

my $homedir = `echo \$HOME`;
my $workdir = getcwd();

my $i=0;
my $j=1;
my $k=1;
my $nlib=150; # number of library spectra close to the global minimum

open (INFILE,"<atmosphericParameters.dat");
my @data = <INFILE>;
my $lines=@data;
close (INFILE);

open (BUFF,">buffer.dat");
for($i=0;$i<$lines;$i=$i+1){
  if($i==($j-1)*($nlib+2)){
      chomp($data[$i]);

#      while($data[$i+$k]=~/-1000.00/){
#	$k=$k+1;
#      }
#      printf BUFF "%30s%s", $data[$i],$data[$i+$k];

      if ($data[$i+1]=~/-1000.00/){
	if ($data[$i+2]=~/-1000.00/){
	  if ($data[$i+3]=~/-1000.00/){
	    printf BUFF "%30s%s", $data[$i],$data[$i+4];
	  }
	  else {printf BUFF "%30s%s", $data[$i],$data[$i+3];}
	}
	else {printf BUFF "%30s%s", $data[$i],$data[$i+2];}
      }
      else {printf BUFF "%30s%s", $data[$i],$data[$i+1];}

      $j=$j+1;
  }
}
close(BUFF);

`awk '{print \$2}' buffer.dat > best-miles.tab`;
