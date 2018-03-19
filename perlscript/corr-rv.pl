#!/usr/bin/perl -w

`ls *.ascii > list`;
`awk '{print \$3}' radialVelocities.dat > vr-only.tab`;
my $i=0;
my $name="";
my @tmp="";
my $size="";
my @vr="";
my $ind=0.0;

open(VR, 'vr-only.tab');
@vr=<VR>;
close(VR);

open (INFILE, 'list');
while(<INFILE>){
chomp;


##########################
# to correct vr from etoile measurements

$name="rv" . $_;
chomp($vr[$i]);
$ind=$vr[$i]*5400.0/299792.458;
`awk '{if(NR==1){printf "%d\\n",\$1}else{printf "%8.2f\\t%9.6f\\n",(\$1*1.0)-($ind*1.0),\$2}}' $_ > $name`;

##########################

print "$_..........ok!\n";
$i=$i+1;
}
close INFILE;

#`rm list`;
