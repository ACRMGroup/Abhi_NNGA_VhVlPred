#! /acrm/usr/local/bin/perl

use strict 'vars';

# ----------- DECLARATION OF GLOBAL VARIABLES ----------

   my $validationPatternFilename="";
   my $neuralNetworkOutputFilename="";
   my $angleComparisonFilename="";

   my $temporaryBatchmanCommandsFilename="";
   my @batchmanCommands=();

   my @actualInterfaceAngles=();
   my @predictedInterfaceAngles=();
   my $interfaceAngleFlag=0;
   my @con=();
   my @values=();
   my $line="";

   my $i=0;

# ---- END OF GLOBAL VARIABLES DECLARATION SECTION -----


# Main code of the program starts here.


if($#ARGV < 2)
{
   print "\nUsage: perl $0 <Arguments>\n";
   print "\nArguments are:\n";
   print "\n1. Path of Neural network output file";
   print "\n2. Path of Validation pattern file";
   print "\n4. Output file to correlate correct and predicted interface angles";
   print "\n\n";
   exit 0;
}

$neuralNetworkOutputFilename=$ARGV[0];
$validationPatternFilename=$ARGV[1];
$angleComparisonFilename=$ARGV[2];

if(! -e $validationPatternFilename)
{
   print "\nFile \"$validationPatternFilename\" doesnt exist.\nExiting program\n\n";
   exit(0);
}

# 1. Read the validation file for correct angles.
#
# 2. Read the neural network output file for predicted angles.
#
# 3. Write actual and corresponding predicted value in seperate columns on the same line
#    so that the two may be correlated.

open(hd,$validationPatternFilename);

foreach $line ( <hd> )
{
   chomp($line);

   if($line =~ /^# PDB/)
   {
      if($#values != -1)
      {
	 push(@actualInterfaceAngles,$values[$#values]);
      }

      @values=();
   }
   elsif($line =~ /^[0-9]/)
   {
      @con=split(/ /,$line);
      push(@values,@con);
   }
}

# Push the last value angle into the array.

push(@actualInterfaceAngles,$values[$#values]);

close(hd);

# Now, perform the same step with the neural network result file.

@con=();
@values=();

open(hd,$neuralNetworkOutputFilename);

foreach $line ( <hd> )
{
   chomp($line);

   if($line =~ /^#/)
   {
      if($#values != -1)
      {
	 push(@predictedInterfaceAngles,$values[$#values]);
      }

      @values=();
   }
   elsif($line =~ /^[0-9]/)
   {
      @con=split(/ /,$line);
      push(@values,@con);
   }
}

close(hd);

# Push the last interface angle into the array.

push(@predictedInterfaceAngles,$values[$#values]);

# Check if the number of actual interfact angles is equal to the number of predicted interface angles
# picked up from the neural network result file. If it isnt, then print an error message.

if($#predictedInterfaceAngles != $#actualInterfaceAngles)
{
   print "\nNumber of interface angles and number of predicted interface angles dont match\n\n";
   exit(0);
}

# Now, write the two arrays into a file.

open(whd,">$angleComparisonFilename");

for($i=0;$i <= $#predictedInterfaceAngles; $i++)
{
   $predictedInterfaceAngles[$i] = sprintf("%3.4f",$predictedInterfaceAngles[$i]);
   $actualInterfaceAngles[$i] = sprintf("%3.4f",$actualInterfaceAngles[$i]);

   print whd "$actualInterfaceAngles[$i]\t$predictedInterfaceAngles[$i]\n";
}

close(whd);

# End of program.
