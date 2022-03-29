#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my %lightHash=();
   my %heavyHash=();
   my $listFilename="";
   my $line="";
   my $pdbCode="";
   my @con=();
   my $floatValueForPosition=-1;
   my @lightFloatingArray=();
   my @heavyFloatingArray=();
   my $position="";
   my $ignore=();
   my $frequency=0;

   my %fval=();

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----

# --------------- SUB - ROUTINES SECTION ---------------

sub get_floating_point_value_for_position
{
   my $position=$_[0];
   my @parts=();
   my $floatValue=0;

   # Remove the first character from the position string -> L35A becomes 35A

   $position=~s/^[A-Z]//;

   # If position doesnt contain an insert code, return the integer value.

   if($position !~ /[A-Z]/)
   {
      # 24

      $floatValue = $position;
   }
   else
   {
      # 35A

      $position=~s/[A-Za-z]//g;

      $floatValue=$position + $fval{$&};
   }

   return $floatValue;
}


# ----------- END OF SUB - ROUTINES SECTION ------------


# Main code of the program starts here.

if($#ARGV < 0)
{
   print "\n Usage: $0 <Input file>\n\n";
   exit(0);
}

$listFilename=$ARGV[0];


# Assign the floating point codes for positions with strings.

$fval{"A"}=0.10;
$fval{"B"}=0.11;
$fval{"C"}=0.12;
$fval{"D"}=0.13;
$fval{"E"}=0.14;
$fval{"F"}=0.15;
$fval{"G"}=0.16;
$fval{"H"}=0.17;
$fval{"I"}=0.18;
$fval{"J"}=0.19;
$fval{"K"}=0.20;
$fval{"L"}=0.21;
$fval{"M"}=0.22;
$fval{"N"}=0.23;
$fval{"O"}=0.24;
$fval{"P"}=0.25;
$fval{"Q"}=0.26;
$fval{"R"}=0.27;
$fval{"S"}=0.28;
$fval{"T"}=0.29;
$fval{"U"}=0.30;
$fval{"V"}=0.31;
$fval{"W"}=0.32;
$fval{"X"}=0.33;
$fval{"Y"}=0.34;
$fval{"Z"}=0.35;

# Open file and read file content into an array.

open(hd,$listFilename);
@con=<hd>;
close(hd);

# Now, process the light and heavy chain positions.

foreach $line (@con)
{
   chomp($line);

   $line=~s/\t/ /;
   $line=~s/\s+/ /;

   if( ($line =~ /^\*L/) || ($line =~ /^L[0-9]/) )
   {
      ($position,$frequency)=split(/ /,$line);
      $position=~s/\*//;
      $floatValueForPosition=&get_floating_point_value_for_position($position);
      $position=sprintf("%5s",$position);
      $frequency=sprintf("%10d",$frequency);
      $lightHash{$floatValueForPosition}="$position $frequency"
   }
   elsif( ($line =~ /^\*H/) || ($line =~ /^H[0-9]/) )
   {
      ($position,$frequency)=split(/ /,$line);
      $position=~s/\*//;
      $floatValueForPosition=&get_floating_point_value_for_position($position);
      $position=sprintf("%5s",$position);
      $frequency=sprintf("%10d",$frequency);
      $heavyHash{$floatValueForPosition}="$position $frequency";
   }
}

# Sort the keys of the two hashes.

@lightFloatingArray=sort { $a <=> $b } (keys %lightHash);
@heavyFloatingArray=sort { $a <=> $b } (keys %heavyHash);

foreach $position (@lightFloatingArray)
{
   print $lightHash{$position},"\n";
}

foreach $position (@heavyFloatingArray)
{
   print $heavyHash{$position},"\n";
}

# End of program.
