# include "intermediate.h"

# define MAXIMUM_NUMBER_OF_INTERFACE_POSITIONS 20

/* ------------------------- GLOBAL VARIABLES DECLARATION SECTION ------------------------- */


/* Variables associated with file names. */

static char interfacePositionsListFilename[100],
	    kabatNumberingFilename[100],
	    residuePropertiesMatrixFilename[100],
	    interfaceAngleFilename[100];

static char pdbCode[8];

static char outputFilename[100];


/* --------------------- END OF GLOBAL VARIABLES DECLARATION SECTION ---------------------- */


/* ---------------------------- FUNCTION DECLARATION SECTION ------------------------------

int read_interface_positions_list(char **interfacePositionsList,char *interfacePositionsListFilename);

KABATLIST *read_kabat_numbered_file(char *kabatNumberingFilename);

BOOL find_residues_at_interface_positions(char **interfacePositionsList,         // Strings of interface positions
                                          int numberOfInterfacePositions,        // Number of elements in the array above
                                          KABATLIST *kabatNumberedListFirstPtr,  // Pointer to the first node of Kabat numbered list
                                          char *residuesAtInterfacePositions);   // Array to which residues must be written

int correlate_interface_positions_and_matrix(FILE *fp,                              // File pointer to file with AA properties
                                             char *residuesAtInterfacePositions,    // List of residues at interface positions
                                             int numberOfInterfacePositions,        // Number of elements in the array above
                                             float **interfaceResiduesMatrix);      // Matrix of Amino acid properties

float extract_interface_angle(FILE *fp,char *pdbCode);

--------------------------- END OF FUNCTION DECLARATION SECTION ------------------------ */


void Usage();

BOOL parse_command_line_parameters(int numberOfParam,char **param);


/* ---------------------- END OF FUNCTION DECLARATION SECTION ----------------------------- */


 
void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nThe arguments are:\n");
   printf("\n1. -ilist <File with list of Interface positions>");
   printf("\n2. -kabnum <File with the Kabat numbered residues>");
   printf("\n3. -mat <File with Matrix representation of the amino acid properties>");
   printf("\n4. -angle <File with interface angles for every antibody>");
   printf("\n5. -pdb <PDB Code of antibody>");
   printf("\n6. -out <Optional parameter - Output file>");
   printf("\n\n");

} /* End of function "Usage". */


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcmp(param[i],"-ilist") )
      {
	 strcpy(interfacePositionsListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-kabnum") )
      {
	 strcpy(kabatNumberingFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-mat") )
      {
	 strcpy(residuePropertiesMatrixFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-angle") )
      {
	 strcpy(interfaceAngleFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-pdb") )
      {
	 strcpy(pdbCode,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-out") )
      {
	 strcpy(outputFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
	 return FALSE;
   }

   return TRUE;

} /* End of function "parse_command_line_parameters". */


int main(int argc,char **argv)
{
   /* Step 1: Read all the command line parameters into program variables */

   if(argc < 11)
   {
      Usage(argv);
      exit(0);
   }

   strcpy(outputFilename,"");

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      exit(0);
   }

   /* Function prototype of create_intermediate_file:

      int create_intermediate_file(char *interfacePositionsListFilename,
                                   char *kabatNumberingFilename,
                                   char *residuePropertiesMatrixFilename,
                                   char *interfaceAngleFilename,
                                   char *pdbCode,
                                   char *outputFilename);
   */

   /* First, read the interface positions list */

   create_intermediate_file(interfacePositionsListFilename,
			    kabatNumberingFilename,
			    residuePropertiesMatrixFilename,
			    interfaceAngleFilename,
			    pdbCode,
			    outputFilename);

   return 1;
}
