#!/usr/bin/perl -w

use Cwd;
#use Env qw(CODE);
#use Switch;

my $homedir = `echo \$HOME`;
my $workdir = getcwd();


my @data;
my $size;
my $i;
my @data2;
my $size2;
my $tmp="";

#`awk '{print \$2}' buffer.dat > best-miles.tab`;

open (STARS, "<best-miles.tab");
@data=<STARS>;
$size=@data;
close (STARS);

open (OUTFILE, ">ref-best-miles.tab");
for($i=0;$i<$size;$i=$i+1){
  chomp($data[$i]);
  $data[$i] =~ s/BD\+/ /; 
  $quit=0;
#  printf OUTFILE "$data[$i]\t";
  open (ALL,"<names_miles91.lis");
  while(<ALL>){
    chomp($_);
    $_ =~ s/BD\+/ /; 
    if($_ =~ m/$data[$i]/ && $quit==0) {
#    if($_ =~ m/HD/ || $_ =~ m/BD-/){$tmp=$_;}
#    else{
#    $_=~ s/^\s+//; #remove leading spaces
#    $tmp="BD+" . $_;}
    print OUTFILE "$_\n"; 
    $quit=1;}
  }
  close(ALL);
}
close(OUTFILE);


`paste best-miles.tab ref-best-miles.tab | awk '{print \$1,\$3}' > ref-best-miles2.tab`;
`awk '{print \$2}' ref-best-miles2.tab > net-only.tab`;
`paste buffer.dat ref-best-miles2.tab | awk '{print \$1,NR,\$9,\$10}' > tmp`;
`wc -l <tmp > header`;
`cat header tmp > templatesForCC.txt`;
`cp -rp templatesForCC.txt ../code/dnm/`;

