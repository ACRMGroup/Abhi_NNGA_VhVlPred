/* This program trains and validates a neural network to predict interface angle from a set
   of interface positions. The set of interface positions is represented by a gene (a series
   of 1s and 0s). A 1 represents that a certain interface position to be included for analysis
   while 0 implies its exclusion.

   The steps involved in the program may be summarised as follows:

   --> Read the required gene from file

   --> Read list of PDB ids.

   --> Read the set of all interface positions into a list.

   --> Read interface angles for all the PDBs into a list.

   --> Read Kabat numbering for the PDBs into a series of linked lists.

   --> Read residue properties into a matrix.

   --> Calculate the correlation coefficient over N folds for the gene and write it into a file
       (Eg. P5_score.txt) along with the gene code and generation in the following format:

       GENE #1: 000101000100101100000000
       SCORE: 0
       GENERATION: 1
*/

# define MAIN 1

# include <stdio.h>
# include <string.h>

# include "grid_ga.h"
# include "ga.h"
/* # include "snns.h" */

/* struct gene_score_generation
   {
      char gene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];
      double score;
      int generation;
   };
*/

struct gene_score_generation geneScoreGen;

/* struct snns_filenames
   {
      // Parameters pertaining to filenames.
   
      char snnsTrainingScriptFilename[MAX_FILENAME_LENGTH],
           snnsValidationScriptFilename[MAX_FILENAME_LENGTH],
           untrainedNetworkFilename[MAX_FILENAME_LENGTH],
           trainingPatternsFilename[MAX_FILENAME_LENGTH],
           validationPatternsFilename[MAX_FILENAME_LENGTH],
           trainedNetworkFilename[MAX_FILENAME_LENGTH],
           predictedValuesOutputFilename[MAX_FILENAME_LENGTH],
           inputOutputComparisonsFilename[MAX_FILENAME_LENGTH],
           coefficientOutputFilename[MAX_FILENAME_LENGTH];
   };
*/

struct snns_filenames snnsFilenames;

char genesFilename[MAX_FILENAME_LENGTH];

int geneNumber=0;

/* struct list_of_pdb
   {
      int numberOfPDB;
      char pdbListFilename[MAX_FILENAME_LENGTH];
      char **pdbList;
   };
*/

struct list_of_pdb listOfPDB;

/* typedef struct interface_positions_and_angles
   {
      // Interface positions
   
      char **interfacePositionsList;
      int numberOfInterfacePositions;
      char interfacePositionsListFilename[MAX_FILENAME_LENGTH];
   
      // Interface angles
   
      double *interfaceAngles;
      char interfaceAnglesFilename[MAX_FILENAME_LENGTH];
   
   }INTERFACE;
*/

INTERFACE interface;

/* typedef struct kabat_list
   {
      char position[8],
           residueOneLetterCode;
   
      struct kabat_list *next;
   
   }KABATLIST;
*/

KABATLIST **kabatNumberingFirstPointers=NULL;

FILE *fp = NULL;

int i=0;

double correlationCoefficientOrError=0;

struct residue_properties residueProperties;

char trainingBoundaryString[20];

/* BOOL writeOutputToFileFlag=FALSE; */

char outputFilename[MAX_FILENAME_LENGTH];

/* Structures used in the program:

   struct gene_score_generation
   {
      char gene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];
      double score;
      int generation;
   };

   struct snns_filenames
   {
      // Parameters pertaining to filenames

      char snnsTrainingScriptFilename[1000],
           snnsValidationScriptFilename[1000],
           untrainedNetworkFilename[1000],
           trainingPatternsFilename[1000],
           validationPatternsFilename[1000],
           trainedNetworkFilename[1000],
           predictedValuesOutputFilename[1000],
           inputOutputComparisonsFilename[1000],
           coefficientOutputFilename[1000];
   };

   struct list_of_pdb
   {
      int numberOfPDB;
      char **pdbList;
      char pdbListFilename[500];
   };

   typedef struct interface_positions_and_angles
   {
      // Interface positions
   
      char **interfacePositionsList;
      int numberOfInterfacePositions;
      char interfacePositionsListFilename[500];
   
      // Interface angles
   
      double *interfaceAngles;
      char interfaceAnglesFilename[500];
   
   }INTERFACE;

   typedef struct kabat_list
   {
      char position[8],
           residueOneLetterCode;
 
      struct kabat_list *next;
 
   }KABATLIST;

*/



/* ----------------------------- FUNCTION DECLARATION SECTION -------------------------- */

void free_memory();

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

/* -------------------------- END OF FUNCTION DECLARATION SECTION ---------------------- */


/* ------------------------------- FUNCTION DEFINITION SECTION ------------------------- */


/* void free_memory(): This function releases memory dynamically allocated during program execution */

void free_memory()
{
   int i=0;

   KABATLIST *p=NULL,
             *prev=NULL;

   /* Structure list_of_pdb

      struct list_of_pdb
      {
         int numberOfPDB;
         char **pdbList;              -> Allocated memory. Number of elements is given by numberOfPDB.
         char pdbListFilename[500];
      }
   */

   for(i=0;i<listOfPDB.numberOfPDB;i++)
   {
      FREE(listOfPDB.pdbList[i]);
   }

   FREE(listOfPDB.pdbList);

   /* Structure interface_positions_and_angles or INTERFACE

         typedef struct interface_positions_and_angles
         {
            // Interface positions
         
            char **interfacePositionsList;            -> Number of elements given by numberOfInterfacePositions
            int numberOfInterfacePositions;
            char interfacePositionsListFilename[500];
         
            // Interface angles
         
            double *interfaceAngles;
            char interfaceAnglesFilename[500];
         
         }INTERFACE;
   */

   for(i=0;i<interface.numberOfInterfacePositions;i++)
   {
      FREE(interface.interfacePositionsList[i]);
   }

   FREE(interface.interfacePositionsList);
   FREE(interface.interfaceAngles);

   /* KABATLIST **kabatNumberingFirstPointers - number of elements is given by listOfPDB.numberOfPDB */

   for(i=0;i<listOfPDB.numberOfPDB;i++)
   {
      p=kabatNumberingFirstPointers[i];

      while(p)
      {
         prev=p;
         p=p->next;
         FREE(prev);
      }
   }

   FREE(kabatNumberingFirstPointers);

} /* End of function "FREE_memory" */


/* void Usage(char **argv):

   This function displays the command line parameters to be input for the program
*/
   
void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nArguments are:\n");
   printf("\n1. -in <File containing genes>");
   printf("\n2. -num <Number of gene to be chosen from file>");
   printf("\n3. -int <File with interface positions>");
   printf("\n4. -pdb <File with list of PDB codes>");
   printf("\n5. -tor <File with interface angles>");
   printf("\n6. -trb <Training boundary string> (Eg: 1-400)");
   printf("\n7. -out <Name of output file>");
   printf("\n\n");

} /* End of function "Usage" */


/* BOOL parse_command_line_parameters(int numberOfParam,char **param):

   This function parses the command line parameters and returns TRUE on successfully parsing
   them. Otherwise, FALSE is returned to the calling function.
*/

BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-in") )
      {
         strcpy(genesFilename,param[i+1]);
         i+=2;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-num") )
      {
         geneNumber=atoi(param[i+1]);
         i+=2;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-int") )
      {
         strcpy(interface.interfacePositionsListFilename,param[i+1]);
         i+=2;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-pdb") )
      {
         strcpy(listOfPDB.pdbListFilename,param[i+1]);
         i+=2;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-tor") )
      {
         strcpy(interface.interfaceAnglesFilename,param[i+1]);
         i+=2;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-trb") )
      {
	 strcpy(trainingBoundaryString,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(outputFilename,param[i+1]);
	 writeOutputToFileFlag = TRUE;
	 i+=2;
	 continue;
      }
      else
      {
         return FALSE;
      }
   }

   return TRUE;

} /* End of function "parse_command_line_parameters" */


/* --------------------------- END OF FUNCTION DEFINITIONS SECTION --------------------- */


int main(int argc,char **argv)
{
   /* Parse command line parameters */

   if(argc < 12)
   {
      Usage(argv);
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      printf("\nERROR PARSING COMMAND LINE PARAMETERS\n");

      Usage(argv);
      return 0;
   }

   /* Check if the environment variable TEMPORARY_SNNS_PATH has been set. If not report error and quit program */

   if(! getenv("TEMPORARY_SNNS_PATH") )
   {
      printf("\nPlease set environment variable \"TEMPORARY_SNNS_PATH\" before executing program\n\n");
      return 0;
   }

   /* Read the required gene from file. The function that does this is:

      read_gene_from_file(char *fileName,int geneNumber,struct gene_score_generation *geneScoreGen);
   */

   if(! read_gene_from_file(genesFilename,geneNumber,&geneScoreGen) )
   {
      printf("\nFile \"%s\" does not exist. Exiting from program\n\n",genesFilename);
      return 0;
   }

   /* Read the list of PDB codes into a list */

   listOfPDB.pdbList=NULL;

   listOfPDB.pdbList=(char **)malloc(MAX_NUMBER_OF_PDB * sizeof(char *));

   listOfPDB.numberOfPDB = read_pdb_list(&listOfPDB);

   /* Read interface positions into a list */

   interface.interfacePositionsList = (char **)malloc(MAX_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char *));

   interface.numberOfInterfacePositions = read_interface_positions_list(interface.interfacePositionsList,
                                                                        interface.interfacePositionsListFilename);

   /* Read the interface angles */

   interface.interfaceAngles=(double *)malloc(listOfPDB.numberOfPDB * sizeof(double));

   fp=fopen(interface.interfaceAnglesFilename,"r");

   if(fp==NULL)
   {
     fprintf(stderr,"File %s not opened for read\n",interface.interfaceAnglesFilename);
     exit(1);
   }
 
   for(i=0;i<listOfPDB.numberOfPDB;i++)
   {
      rewind(fp);

      interface.interfaceAngles[i]=extract_interface_angle(fp,listOfPDB.pdbList[i]);
   }

   fclose(fp);

   /* Create an array of linked lists for the Kabat numbering of all antibodies.
      This is done by calling the function:

      BOOL create_kabat_numbering_lists(struct list_of_pdb *listOfPDB,
                                        int numberOfPDB,
                                        KABATLIST **kabatNumberingFirstPointers);
   */

   kabatNumberingFirstPointers = (KABATLIST **)malloc(MAX_NUMBER_OF_PDB * sizeof(KABATLIST *));

   create_kabat_numbering_lists(&listOfPDB,listOfPDB.numberOfPDB,kabatNumberingFirstPointers);

   /* Assign filenames to all files in "struct snns_filenames snnsFilenames".
      The format of this is as follows:

      struct snns_filenames
      {
         char snnsTrainingScriptFilename[1000],
              snnsValidationScriptFilename[1000],
              untrainedNetworkFilename[1000],
              trainingPatternsFilename[1000],
              validationPatternsFilename[1000],
              trainedNetworkFilename[1000],
              predictedValuesOutputFilename[1000],
              inputOutputComparisonsFilename[1000],
              coefficientOutputFilename[1000];
      };

      Prototype of the function that assigns the file names:

      void assign_snns_filenames(struct snns_filenames *snnsFilenames,
                                 int currentProcessID,
                                 BOOL createUntrainedNetworkFileFlag); // Flag to set untrained network file name in function.
   */

   assign_snns_filenames(&snnsFilenames,
                         getpid(),
                         TRUE);

   /* If the output is to be written to a file (indicated by the flag writeOutputToFileFlag), then
      copy the string in outputFilename to snnsFilenames.coefficientOutputFilename.
   */

   if(writeOutputToFileFlag)
   {
      strcpy(snnsFilenames.overallPredictionsFilename, outputFilename);
   }

   /* Read the residue properties into "struct residue_properties".

      struct residue_properties
      {
         char *residues;
         int numberOfResiduesInMatrix;
         double **residuePropertiesMatrix;
         char residuePropertiesMatrixFilename[500];
      };

      The function that does this is:

      void read_residue_properties_matrix_into_struct(struct residue_properties *residueProperties);
   */

   sprintf(residueProperties.residuePropertiesMatrixFilename,"%s/ABNUM/residue_matrix.mat",getenv("ABHIDATADIR"));

   if(! access(residueProperties.residuePropertiesMatrixFilename, R_OK) )
   {
      fprintf(stderr, "\nUnable to open file \"%s\" in read mode.\n\n", residueProperties.residuePropertiesMatrixFilename);
      return(0);
   }

   read_residue_properties_matrix_into_struct(&residueProperties);

   /* Now, calculate coefficient over the folds. This is done using the following function:

      double calculate_coefficient_over_single_fold(char *gene,
                                                    char *fold,
                                                    INTERFACE *interface,
                                                    struct residue_properties *residueProperties,
                                                    struct snns_filenames *snnsFilenames,
                                                    KABATLIST **kabatNumberingFirstPointers,
                                                    struct list_of_pdb *listOfPDB,
                                                    char *desc);
   */

   correlationCoefficientOrError = calculate_score_over_single_fold(geneScoreGen.gene,
                                                                    trainingBoundaryString,
                                                                    &interface,
                                                                    &residueProperties,
                                                                    &snnsFilenames,
                                                                    kabatNumberingFirstPointers,
                                                                    &listOfPDB,
                                                                    "");

   /* Print correlation coefficient or error depending on the status of the flag isCorrelationCoefficient */

   if(isCorrelationCoefficient)
   {
      printf("\nCorrelation coefficient: %4.4f\n",correlationCoefficientOrError);
   }
   else
   {
      printf("\nError in prediction of interface angle fraction: %4.4f\n",correlationCoefficientOrError);
   }

   /* Return 1 and exit from program */

   free_memory();

   return 1;

} /* End of program */