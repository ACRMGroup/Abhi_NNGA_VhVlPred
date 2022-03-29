# include "ga.h"
# define NUMBER_OF_FOLDS 5

/* Global variables

   struct gene
   {
      char **parentGenes;
      char **childGenes;
      char *bestGene;
   
      int lengthOfGene;
      double *parentScores;
      double *childScores;
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

   struct list_of_pdb
   {
      int numberOfPDB;
      char **pdbList;
      char pdbListFilename[500];
   };

   struct residue_properties
   {
      char *residues;
      int numberOfResiduesInMatrix;
      double **residuePropertiesMatrix;
      char residuePropertiesMatrixFilename[500];
   };
*/

struct gene population;
INTERFACE interface;
struct list_of_pdb listOfPDB;
struct residue_properties residueProperties;

 
/* OTHER GENERAL VARIABLES */

char **folds=NULL;

static int randomNumberSeed=0; /* Seed for the random number generator */

double overallScores[2 * MAX_POPULATION]; /* Array used to store scores of the parent and child */
					  /* genes while selecting the best amongst them */

KABATLIST **kabatNumberingFirstPointers=NULL;


/* ------------------------- FUNCTION DECLARATION SECTION ------------------------- */

BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-pdb") )
      {
	 strcpy(listOfPDB.pdbListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-pos") )
      {
	 strcpy(interface.interfacePositionsListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-ang") )
      {
	 strcpy(interface.interfaceAnglesFilename,param[i+1]);
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


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n\nArguments are:\n",argv[0]);
   printf("\n1. -pdb <File with list of PDB codes>");
   printf("\n2. -pos <File with list of kabat numbered interface positions>");
   printf("\n3. -ang <File with list of interface angle listed according to PDB code");
   printf("\n\n");

} /* End of function "Usage" */


/* ---------------------- END OF FUNCTION DEFINITION SECTION ---------------------- */

int main(int argc,char **argv)
{
   FILE *fp=NULL;

   int i=0;

   double bestScore=1;

   struct snns_filenames snnsFilenames;

   /* Initialise the random number seed */

   randomNumberSeed=getpid();

   /* Check if sufficient number of arguments have been provided */

   if(argc < 6)
   {
      Usage(argv);
      return 0;
   }

   /* Parse the command line parameters */

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   /* Check if the environment variable TEMPORARY_SNNS_PATH has been set. If not report error and quit program */

   if(! getenv("TEMPORARY_SNNS_PATH") )
   {
      printf("\nPlease set environment variable \"TEMPORARY_SNNS_PATH\" before executing program\n\n");
      return 0;
   }

   /* Read the list of PDB codes into a list -> COMPLETED */

   listOfPDB.pdbList=(char **)malloc(1000 * sizeof(char *));

   listOfPDB.numberOfPDB = read_pdb_list(&listOfPDB);

   /* Create the folds for training and validation of the network -> COMPLETED */

   folds=(char **)malloc(NUMBER_OF_FOLDS * sizeof(char *));

   create_fold_ranges(folds,listOfPDB.numberOfPDB,NUMBER_OF_FOLDS);

   /* Read interface positions into a list -> Available in the library */

   interface.interfacePositionsList = (char **)malloc(MAX_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char *));

   interface.numberOfInterfacePositions = read_interface_positions_list(interface.interfacePositionsList,
									interface.interfacePositionsListFilename);

   /* Read the interface angles */

   interface.interfaceAngles=(double *)malloc(listOfPDB.numberOfPDB * sizeof(double));

   fp=fopen(interface.interfaceAnglesFilename,"r");

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


   kabatNumberingFirstPointers=(KABATLIST **)malloc(MAX_NUMBER_OF_PDB * sizeof(KABATLIST *));

   create_kabat_numbering_lists(&listOfPDB,listOfPDB.numberOfPDB,kabatNumberingFirstPointers);

   /* Allocate space for the parent and child genes */

   population.parentGenes=(char **)malloc(POPULATION * sizeof(char *));
   population.childGenes=(char **)malloc(POPULATION * sizeof(char *));

   for(i=0;i<POPULATION;i++)
   {
      population.parentGenes[i]=(char *)malloc((numberOfInterfacePositions + 1) * sizeof(char));
      population.childGenes[i]=(char *)malloc((numberOfInterfacePositions + 1) * sizeof(char));
   }

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

   /* Create the initial set of Parent genes -> COMPLETED */

   initialise_population(population.parentGenes,POPULATION,interface.numberOfInterfacePositions);

   /* While overall correlation coefficient of the gene set is less than 0.75, do the following steps */

   bestScore=0;

   population.parentScores=(double *)malloc(POPULATION * sizeof(double));
   population.childScores=(double *)malloc(POPULATION * sizeof(double));

   population.bestGene=(char *)malloc( (numberOfInterfacePositions + 1) * sizeof(char));

   /* Calculate the correlation coefficient for the parent genes using the function:

      double calculate_coefficient_over_folds(char *gene,
                                              char **folds,
                                              int numberOfFolds,
                                              INTERFACE *interface,
                                              struct residue_properties *residueProperties,
                                              struct snns_filenames *snnsFilenames,
                                              KABATLIST **kabatNumberingFirstPointers,
                                              struct list_of_pdb *listOfPDB);
   */

   i=0;

   while(i < POPULATION)
   {
      population.parentScores[i] = calculate_coefficient_over_folds(population.parentGenes[i],
								    folds,
								    NUMBER_OF_FOLDS,
								    &interface,
								    &residueProperties,
								    &snnsFilenames,
								    kabatNumberingFirstPointers,
								    &listOfPDB);
      i++;
   }

   while(bestScore < 0.75)
   {
      /* Create a child gene set with POPULATION genes. This is done using the function:

	 void create_genes(char **parentGenes,
	                   double *parentScores,
	                   char **childGenes,
	                   int population,
	                   int lengthOfGene,
	                   int numberOfONPositions);
      */

      create_genes(population.parentGenes,
		   population.parentScores,
		   population.childGenes,
		   POPULATION,
		   interface.numberOfInterfacePositions,
		   NUMBER_OF_ON_POSITIONS);

      /* Calculate scores for the child genes */

      i=0;

      while(i < POPULATION)
      {
	 childScores[i]=calculate_coefficient_over_folds(childGenes[i],
							 folds,
							 interfacePositionsList,
							 numberOfInterfacePositions,
							 &snnsFilenames);
	 i++;
      }

      /* Now choose a new gene set with POPULATION genes which have the best scores */

      bestScore=select_best_genes(parentGenes,parentScores,
				  childGenes,childScores,
				  POPULATION,numberOfInterfacePositions,
				  bestGene);

      printf("\nBest gene: %s\nBest score: %f",bestGene,bestScore);
      fflush(stdout);

   } /* End of while loop */


   /* Remove files and release memory */

   remove_snns_files(&snnsFilenames);

   free_array_2D_char(parentGenes,POPULATION);
   free_array_2D_char(childGenes,POPULATION);

   free(bestGene);

   free(parentScores);
   free(childScores);

   free_array_2D_char(interfacePositionsList,numberOfInterfacePositions);

   free(interfaceAngles);

   free_array_2D_char(folds,NUMBER_OF_FOLDS);

   free(residues);

   free_array_2D_double(residuePropertiesMatrix,numberOfResiduesInMatrix);

   for(i=0;i<numberOfPDB;i++)    /* Release the Kabat numbered lists */
   {
      p=kabatNumberingFirstPointers[i];

      while(p)
      {
	 prev=p;
	 p=p->next;
	 free(prev);
      }
   }

   /* Return 1 to the compiler */

   return 1;

} /* End of program */
