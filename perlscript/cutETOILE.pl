#!/usr/bin/perl -w

`ls *.ascii > list`;
my $i=0;
my $name="";
my $size="";
my $ind=0.0;
my $id="";

open (INFILE, 'list');
while(<INFILE>){
chomp;

##########################
# to clean and put the number of lines at the header
$id=substr($_,5,6);
$name="Ter9-" . $id . ".ascii";
`awk '{if(\$2>0.01 && \$1<6.000*1E+3 && \$1>4.800*1E+3){printf\"%7.1f\\t  %9.4f\\n\",\$1,\$2}}' $_ > tmp`;

`wc -l tmp | awk '{print \$1}' > header`;
`cat header tmp > $name`;

`rm header`;
`rm tmp`;

##########################

print "$_..........ok!\n";
$i=$i+1;
}
close INFILE;

`rm list`;

#Terzan9specid000000245jd2457544p7108f014.fits
