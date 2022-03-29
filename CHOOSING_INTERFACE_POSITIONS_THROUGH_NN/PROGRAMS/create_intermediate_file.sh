# Step 1: Check environment variable ABINTERFACE is set to appropriate directory.

if [ ! `echo $ABINTERFACE` ]
then
   echo
   echo "Please set environment variable \$ABINTERFACE to the appropriate directory"
   echo
   exit 0
fi

# Step 2: Check if file "Fv_Fab_pdb.lst" exists in the directory

if [ ! `ls $ABINTERFACE/Fv_Fab_pdb.lst` ]
then
   echo
   echo "File \"Fv_Fab_pdb.lst\" missing from directory $ABINTERFACE"
   echo
   exit 0
fi

# Step 3: Check number of command line parameters. Display usage if needed.

if [ $# -lt 2 ]
then
   echo
   echo "Usage: sh $0 <Arguments>"
   echo
   echo "Arguments are:"
   echo
   echo "1. File containing list of interface positions"
   echo "2. Name of output file (Intermediate file)"
   echo
   exit 0
fi

fileWithInterfacePositions=$1
outputFilename=$2

residueMatrixFilename=$ABHIDATADIR/ABNUM/residue_matrix.mat

# Step 4: Check if file with interface positions exists. If it doesnt, exit from program.

if [ ! -e $fileWithInterfacePositions ]
then
   echo
   echo "File \"$fileWithInterfacePositions\" missing. Exiting program."
   echo
   exit 0
fi

if [ ! -e $residueMatrixFilename ]
then
   echo
   echo "File \"$residueMatrixFilename\" does not exist"
   echo
   exit 0
fi

# Step 5: Execute the program "create_intermediate_file" to link all the data files and produce an output file.

rm -f $outputFilename

for pdbCode in `cat $ABINTERFACE/Fv_Fab_pdb.lst`
do
   echo $pdbCode

   command="./create_intermediate_file.exe -ilist $fileWithInterfacePositions -kabnum $ABINTERFACE/NUMBERED_FILES/$pdbCode.pro -mat $residueMatrixFilename -angle $ABINTERFACE/INTERFACE_ANGLE/torsion_angles.txt -pdb $pdbCode -out $outputFilename"

   echo $command
   read

done
