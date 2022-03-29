# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>

# include "random.h"
# include "intermediate.h"
# include "snns.h"
# include "arrays.h"

# define POPULATION 10
# define MAX_POPULATION 10
# define MAX_NUMBER_OF_INTERFACE_POSITIONS 200
# define NUMBER_OF_FOLDS 5
# define MAX_NUMBER_OF_PDB 1000
# define TEMPORARY_FILES_DIRECTORY "/home/bsm2/abhi/SNNS/tmp"
# define PC 0.7
# define NUMBER_OF_HIDDEN_NODES 10
# define NUMBER_OF_ON_POSITIONS 20

/* GENES */

static char **parentGenes=NULL,
	    **childGenes=NULL,
	    *bestGene;

static double *parentScores=NULL,
	      *childScores=NULL;

/* INTERFACE POSITIONS AND ANGLES */

static char **interfacePositionsList=NULL; /* A 2D array where all interface positions are listed.
					      This in turn reflects the order of positions in the genes */

static char interfacePositionsListFilename[200];

static int numberOfInterfacePositions=0;

static char interfaceAnglesFilename[200];

static double *interfaceAngles=NULL;

/* KABAT NUMBERING OF PDB FILES */

static KABATLIST **kabatNumberingFirstPointers=NULL;

static KABATLIST *p=NULL,     /* Pointers used to release memory allocated to the */
		 *prev=NULL;  /* Kabat numbering linked lists for the PDB files */

static char **folds=NULL;

/* PDB */

extern int numberOfPDB;

static char **pdbList=NULL;

static char pdbListFilename[100];

/* Declaration of variables for the residue properties */

static char *residues=NULL;

static int numberOfResiduesInMatrix=0;

static double **residuePropertiesMatrix=NULL;

static char residuePropertiesMatrixFilename[200];

/* OTHER GENERAL VARIABLES */

static int randomNumberSeed=0; /* Seed for the random number generator */

double overallScores[2 * MAX_POPULATION]; /* Array used to store scores of the parent and child */
					  /* genes while selecting the best amongst them */


/* ------------------------- FUNCTION DECLARATION SECTION ------------------------- */

int compare(const void *a, const void *b);

int number_of_occurrences(char *string,char ch);

void initialise_population(char **parentGenes,int population,int lengthOfGene,int numberOfONPositions);

int read_pdb_list(char *pdbListFilename,char **pdbList);

void create_fold_ranges(char **folds,int numberOfPDB,int numberOfFolds);

BOOL create_kabat_numbering_lists(char **pdbList,int numberOfPDB,KABATLIST **kabatNumberingFirstPointers);

double calculate_coefficient_over_folds(char *gene,
                                        char **folds,
                                        char **interfacePositionsList,
                                        int numberOfInterfacePositions,
					struct snns_filenames *snnsFilenames);

double calculate_coefficient_over_single_fold(char *gene,
                                              char *fold,
                                              char **interfacePositionsList,
                                              int numberOfInterfacePositions,
                                              struct snns_filenames *snnsFilenames);

void create_genes(char **parentGenes,
                  double *parentScores,
                  char **childGenes,
                  int population,
                  int lengthOfGene,
		  int numberOfONPositions);

double select_best_genes(char **parentGenes,double *parentScores,
                         char **childGenes,double *childScores,
                         int population,int numberOfInterfacePositions,
			 char *bestGene);

BOOL parse_command_line_parameters(int param,char **numberOfParam);

void Usage(char **argv);

/* --------------------- END OF FUNCTION DECLARATION SECTION ---------------------- */



/* -------------------------- FUNCTION DEFINITION SECTION ------------------------- */


int compare(const void *a, const void *b)
{
   if(overallScores[*(int *)a] < overallScores[*(int *)b])
      return(1);
   else if(overallScores[*(int *)a] == overallScores[*(int *)b])
      return(0);
   else
      return(-1);

} /* End of function "compare" */


int number_of_occurrences(char *string,char ch)
{
   int count=0,
       i=0;

   for(i=0;i<strlen(string);i++)
   {
      if(string[i] == ch)
      {
	 count++;
      }
   }

   return count;

} /* End of function "number_of_occurrences" */



void initialise_population(char **parentGenes,	    /* Parent genes */
			   int population,	    /* Number of genes */
			   int lengthOfGene,	    /* Length of gene */
			   int numberOfONPositions) /* Number of positions to switch on (1) */
{
   int i=0,j=0;
   double randomFloatingPointNumber;

   int positionInGene=0;

   for(i=0;i<population;i++)
   {
      for(j=0;j<lengthOfGene;j++)
      {
	 randomFloatingPointNumber=get_random_number_double(randomNumberSeed++);

	 if(randomFloatingPointNumber > PC)
	 {
	    parentGenes[i][j]='1';
	 }
	 else
	 {
	    parentGenes[i][j]='0';
	 }
      }

      parentGenes[i][lengthOfGene]='\0';
   }

   for(i=0;i < population;i++)
   {
      if( number_of_occurrences(parentGenes[i],'1') <= numberOfONPositions )
      {
	 continue;
      }

      /* We have to reduce the number of 1's in the gene */

      while( number_of_occurrences(parentGenes[i],'1') > numberOfONPositions )
      {
	 positionInGene=0;

	 while(parentGenes[i][positionInGene] != '1')
	 {
	    positionInGene=(int)(get_random_number_double(randomNumberSeed++) * lengthOfGene);
	 }

	 parentGenes[i][positionInGene]='0';
      }
   }

} /* End of function "initialise_population" */


int read_pdb_list(char *pdbListFilename,char **pdbList)
{
   FILE *fp=fopen(pdbListFilename,"r");

   char line[100];

   int i=0;

   if(! fp)
   {
      return -1; /* Return -1 if file does not exist */
   }

   while( fgets(line,100,fp) )
   {
      line[strlen(line)-1]='\0';

      pdbList[i]=(char *)malloc((strlen(line)+1) * sizeof(char));

      strcpy(pdbList[i++],line);
   }

   fclose(fp);

   return i;

} /* End of function "read_pdb_list" */


void create_fold_ranges(char **folds,int numberOfPDB,int numberOfFolds)
{
   int pdbPerFold=numberOfPDB/numberOfFolds;

   int i=0;

   int numberOfPDBPerTrainingSet=(numberOfFolds-1)*pdbPerFold;

   int lowerBoundary=0,upperBoundary=0;

   for(i=1;i<=numberOfFolds;i++)
   {
      lowerBoundary = ( (i-1) * pdbPerFold ) + 1;
      upperBoundary = lowerBoundary + numberOfPDBPerTrainingSet;

      folds[i-1]=(char *)malloc(16 * sizeof(char *));

      if(upperBoundary > numberOfPDB)
      {
	 upperBoundary-=numberOfPDB;
      }

      sprintf(folds[i-1],"%d-%d",lowerBoundary,upperBoundary);
   }

} /* End of function "create_fold_ranges" */


BOOL create_kabat_numbering_lists(char **pdbList,int numberOfPDB,KABATLIST **kabatNumberingFirstPointers)
{
   int i=0;

   char kabatNumberingFilename[100];

   char *ABINTERFACE;

   ABINTERFACE=getenv("ABINTERFACE");

   if(! ABINTERFACE)
   {
      return FALSE;
   }

   for(i=0;i<numberOfPDB;i++)
   {
      sprintf(kabatNumberingFilename,"%s/NUMBERED_FILES/%s.pro",ABINTERFACE,pdbList[i]);

      kabatNumberingFirstPointers[i]=read_kabat_numbered_file(kabatNumberingFilename);
   }

   return TRUE;

} /* End of function "create_kabat_numbering_lists". */


double calculate_coefficient_over_folds(char *gene,
                                        char **folds,
                                        char **interfacePositionsList,
                                        int numberOfInterfacePositions,
					struct snns_filenames *snnsFilenames)
{
   /* The steps involved in calculating the coefficient over several folds is as follows:

      1. From the set of 1's in the gene, create a list of interface positions.		|
											|
      2. Create an intermediate file for the set of interface positions.		|   These steps are
											\   implemented in
      3. Create the training and validation pattern file.				 \  the function
											  \ "calculate_coefficient_over_single_fold"
      4. Create an untrained network with the requisite number of nodes.		  /
											 /
      5. SNNS part: Write a script that loads network, trains it on the pattern file,	|
		    evaluates it on the validation file, and returns a coefficient.	|
											|
      6. Do the above steps for every fold.

      7. Calculate the coefficient over folds as being the average of all the coefficients.
   */

   int i=0;

   double pearsonsCoefficients[NUMBER_OF_FOLDS],
	  total=0;

   /* Evaluate the correlation coefficient over every fold and get the average.
      Prototype of the function that does this is:

      double calculate_coefficient_over_single_fold(char *gene,
						    char *fold,
						    char **interfacePositionsList,
						    int numberOfInterfacePositions,
						    struct snns_filenames *snnsFilenames);
   */

   total=0;

   for(i=0;i<NUMBER_OF_FOLDS;i++)
   {
      pearsonsCoefficients[i] = calculate_coefficient_over_single_fold(gene,
								       folds[i],
								       interfacePositionsList,
								       numberOfInterfacePositions,
								       snnsFilenames);

      total+=pearsonsCoefficients[i];
   }

   /* Return the average coefficient over all the folds */

   return total/NUMBER_OF_FOLDS;

} /* End of function "calculate_coefficient_over_folds". */



double calculate_coefficient_over_single_fold(char *gene,
					      char *fold,
					      char **interfacePositionsList,
					      int numberOfInterfacePositions,
					      struct snns_filenames *snnsFilenames)
{
   /* Steps involved in this function

      1. From the set of 1's in the gene, create a list of interface positions.

      2. Create an intermediate file for the set of interface positions.

      3. Create the training and validation pattern file.

      4. Create an untrained network with the requisite number of nodes.

      5. SNNS part: Write a script that loads network, trains it on the pattern file,
                    evaluates it on the validation file, and returns a coefficient.
   */


   /* ############# DECLARATION FOR LOCAL VARIABLES ############# */

   /* INTERFACE POSITIONS AND INTERFACE RESIDUES */

   char **currentInterfacePositions=NULL;

   int currentNumberOfInterfacePositions=0;

   char *residuesAtInterfacePositions=NULL;

   double **interfaceResiduesMatrix=NULL;

   /* INTERMEDIATE FILE */

   char intermediateFilename[200];

   /* TRAINING AND VALIDATION PATTERNS */

   struct pattern_variables patternVariables;

   struct snns_pattern_entry *trainingStartPointer=NULL,
			     *validationStartPointer=NULL;

   struct snns_pattern_entry *p=NULL,	  /* Pointers to delete nodes in the training and */
			     *prev=NULL;  /* validation pattern linked lists */

   /* NEURAL NETWORK */

   struct snns_training_parameters snnsTrainingParameters;

   struct untrained_network untrainedNetwork;

   double pearsonsCoefficient=0;

   /* GENERAL VARIABLES */

   FILE *wfp=NULL;

   BOOL isCorrelationCoefficient=FALSE;

   int i=0,j=0;

   int processID=0;

   /* #################### END OF LOCAL VARIABLES ############### */


   /* Create a list of current interface positions by looking at the gene.
      All interface positions that correspond to a 1 in the gene must be included
      for current examination.
   */

   currentInterfacePositions=(char **)malloc(MAX_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char *));

   for(i=0;i<numberOfInterfacePositions;i++)
   {
      if(gene[i] == '1')
      {
	 currentInterfacePositions[j]=(char *)malloc( (strlen(interfacePositionsList[i]) + 1) * sizeof(char *));
	 strcpy(currentInterfacePositions[j++],interfacePositionsList[i]);
      }
   }

   currentNumberOfInterfacePositions=j;

   /* Create the intermediate file for the given list of interface positions.
      This is done by appending to the intermediate file, the inputs and output
      for every PDB entry.

      The prototype of the function to create an intermediate file is as follows:

      void do_create_intermediate_file(char **interfacePositionsList,
                                       int numberOfInterfacePositions,
                                       char *residuesAtInterfacePositions,
                                       double **interfaceResiduesMatrix,
                                       double interfaceAngle,
                                       char *pdbCode,
                                       FILE *wfp)
   */

   processID=getpid();
   sprintf(intermediateFilename,"%s/%d_intermediate.dat",TEMPORARY_FILES_DIRECTORY,processID);
   wfp=fopen(intermediateFilename,"w");

   sprintf(residuePropertiesMatrixFilename,"%s/ABNUM/residue_matrix.mat",getenv("ABHIDATADIR"));

   /* Read the matrix represnting the amino acids.

      int read_residue_properties_matrix(char *matrixFilename,                // Name of matrix file
                                         int numberOfColumns,                 // Number of columns representing the AA properties
                                         char *residues,                      // Array that will contain order of residues read
                                         double **residuePropertiesMatrix)    // The actual residue property matrix

      This function returns the number of residues that were successfully read from the matrix file.
   */

   residuePropertiesMatrix=(double **)malloc(25 * sizeof(double *));
   residues=(char *)malloc(25 * sizeof(char));

   numberOfResiduesInMatrix = read_residue_properties_matrix(residuePropertiesMatrixFilename,
							     NUMBER_OF_RESIDUE_PROPERTIES_IN_MATRIX,
							     residues,
							     residuePropertiesMatrix);

   if(numberOfResiduesInMatrix < 0)
   {
      return -1000;
   }

   /* Allocate memory for the residues at the interface positions and the matrix that represents them */

   residuesAtInterfacePositions=(char *)malloc( (currentNumberOfInterfacePositions + 1) * sizeof(char));
   interfaceResiduesMatrix=(double **)malloc(currentNumberOfInterfacePositions * sizeof(double *));

   for(i=0;i<numberOfPDB;i++)
   {
      /* Find the residues at the interface positions for the current antibody structure.
	 The function prototype is:

	 BOOL find_residues_at_interface_positions(char **interfacePositionsList,
	                                           int numberOfInterfacePositions,
	                                           KABATLIST *kabatNumberedListFirstPtr,
	                                           char *residuesAtInterfacePositions);
      */

      find_residues_at_interface_positions(currentInterfacePositions,
					   currentNumberOfInterfacePositions,
					   kabatNumberingFirstPointers[i],
					   residuesAtInterfacePositions);

      /* Generate the matrix to represent the residues at the interface positions.
	 Function prototype:

	 int do_correlate_interface_positions_and_matrix(char *residues,
							 double **residuePropertiesMatrix,
							 int numberOfResiduesInMatrix,
							 char *residuesAtInterfacePositions,
							 int numberOfInterfacePositions,
							 double **interfaceResiduesMatrix)
      */

      do_correlate_interface_positions_and_matrix(residues,
						  residuePropertiesMatrix,
						  numberOfResiduesInMatrix,
						  residuesAtInterfacePositions,
						  currentNumberOfInterfacePositions,
						  interfaceResiduesMatrix);

      /* Interface angle is present in the array interfaceAngles

         Now, write intermediate file:

         void do_create_intermediate_file(char **interfacePositionsList,
					  int numberOfInterfacePositions,
					  char *residuesAtInterfacePositions,
					  double **interfaceResiduesMatrix,
					  double interfaceAngle,
					  char *pdbCode,
					  FILE *wfp)
      */

      do_create_intermediate_file(currentInterfacePositions,
				  currentNumberOfInterfacePositions,
				  residuesAtInterfacePositions,
				  interfaceResiduesMatrix,
				  interfaceAngles[i],
				  pdbList[i],
				  wfp);

   } /* End of for(i=0;i<numberOfPDB;i++) */

   fclose(wfp); /* Writing the intermediate file is complete */

   /* Create the training set and validation set lists by reading the intermediate file.

      Prototype of the function for doing this is as follows:

      int read_intermediate_file(char *intermediateFilename,
				 struct pattern_variables *patternVariables,
				 struct snns_pattern_entry **trainingStartPointer,
				 struct snns_pattern_entry **validationStartPointer)

      The definition for "struct pattern_variables" is:

      struct pattern_variables
      {
         char trainingSetBoundaryString[8];
      
         int numberOfTrainingPatterns,
             numberOfValidationPatterns;
      };
   */

   strcpy(patternVariables.trainingSetBoundaryString,fold);

   read_intermediate_file(intermediateFilename,
			  &patternVariables,
			  &trainingStartPointer,
			  &validationStartPointer);

   /* Now, write the training set and validation set into two seperate files.
      Prototype of function to be used:

      void write_patterns_into_file(FILE *wfp,
                                    struct snns_pattern_entry *startPointer,
                                    int numberOfPatterns,
                                    int numberOfInputUnits,
                                    int numberOfOutputUnits);
   */

   wfp=fopen(snnsFilenames->trainingPatternsFilename,"w"); /* Write training patterns file */

   write_patterns_into_file(wfp,
			    trainingStartPointer,
			    patternVariables.numberOfTrainingPatterns,
			    currentNumberOfInterfacePositions * NUMBER_OF_RESIDUE_PARAMETERS,
			    1);

   fclose(wfp);

   wfp=fopen(snnsFilenames->validationPatternsFilename,"w"); /* Write validation patterns file */

   write_patterns_into_file(wfp,
			    validationStartPointer,
			    patternVariables.numberOfValidationPatterns,
			    currentNumberOfInterfacePositions * NUMBER_OF_RESIDUE_PARAMETERS,
			    1);

   fclose(wfp);

   /* Create an untrained network with the appropriate number of input and output units.
      An untrained network is represented by the following data structure.

      struct untrained_network
      {
         int numberOfInputColumns,numberOfInputRows;
         int numberOfHiddenColumns,numberOfHiddenRows;
         int numberOfOutputColumns,numberOfOutputRows;
      };

      The members of this structure are initialized using the function:

      void initialise_untrained_network_parameters(struct untrained_network *untrainedNetwork,
                                                   int numberOfInputRows,int numberOfInputColumns,
                                                   int numberOfHiddenRows,int numberOfHiddenColumns,
                                                   int numberOfOutputRows,int numberOfOutputColumns);

      The network is created using the following function:

      int create_untrained_network(char *untrainedNetworkFilename,struct untrained_network *untrainedNetwork);
   */

   initialise_untrained_network_parameters(&untrainedNetwork,
					   currentNumberOfInterfacePositions,NUMBER_OF_RESIDUE_PARAMETERS,
					   NUMBER_OF_HIDDEN_NODES,1,
					   1,1);

   create_untrained_network(snnsFilenames->untrainedNetworkFilename,&untrainedNetwork);


   /* Write scripts to train and validate the network. Prototype of the function that does this:

      int write_snns_commands(struct snns_filenames *snnsFilenames,
                              struct snns_training_parameters *snnsTrainingParameters);

      Before writing the commands, the training parameters of the neural network have to be initialized.
      The function that does the training parameter initialization is declared the following way:

      BOOL initialise_snns_training_parameters(char *configFilename,struct snns_training_parameters *snnsTrainingParameters);
   */

   if(! initialise_snns_training_parameters(NULL,&snnsTrainingParameters) )
   {
      return -2;
   }

   if(! write_snns_commands(snnsFilenames,&snnsTrainingParameters) )
   {
      return -3;
   }

   /* Execute script to train and validate the neural network. Do the following steps:

      1. Train the network and have it predict angle values for the entries in the validation set.
      2. Get the Pearsons coefficient for the predicted and the actual angle values.

      The function for doing this is declared as follows:

      double execute_script_and_calculate_correlation_coefficient(struct snns_pattern_entry *validationFirst,
                                                                  struct snns_filenames *snnsFilenames,
                                                                  struct pattern_variables *patternVariables,
                                                                  BOOL *isCorrelationCoefficient,
                                                                  BOOL writeResultIntoFileFlag);
   */

   pearsonsCoefficient=execute_script_and_calculate_correlation_coefficient(validationStartPointer,
									    snnsFilenames,
									    &patternVariables,
									    &isCorrelationCoefficient,
									    FALSE);

   /* Release memory */

   free_array_2D_char(currentInterfacePositions,currentNumberOfInterfacePositions);
   free(residuesAtInterfacePositions);
   free_array_2D_double(interfaceResiduesMatrix,currentNumberOfInterfacePositions);

   p=trainingStartPointer;  /* Free the training patterns list */

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   p=validationStartPointer; /* Free the validation patterns list */

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   /* Now return the correlation coefficient to the calling routine */

   return pearsonsCoefficient;

} /* End of function "calculate_coefficient_over_single_fold" */


void create_genes(char **parentGenes,
		  double *parentScores,
		  char **childGenes,
		  int population,
		  int lengthOfGene,
		  int numberOfONPositions)
{
   int i=0;

   int augmentedParentScores[MAX_POPULATION],
       totalScore=0;

   int g1=0,
       g2=0;

   BOOL g1Complete=FALSE,
	g2Complete=FALSE;

   int numberOfChildGenes=0;

   int total=0,
       positionInGene=0;

   char *p=NULL;

   double randomNumber=0;

   /* Scale the scores of the parent genes and total them */

   for(i=0;i<population;i++)
   {
      augmentedParentScores[i]=(int)(parentScores[i] * 1000);
      totalScore+=augmentedParentScores[i];
   }

   /* Perform the following steps till there are population genes in the child set */

   numberOfChildGenes=0;

   while(numberOfChildGenes < population)
   {
      /* Choose two genes by picking random numbers */

      g1 = get_random_number_double(randomNumberSeed++) * totalScore;
      g2 = get_random_number_double(randomNumberSeed++) * totalScore;

      /* Check where in the scale the two numbers g1 and g2 belong.
         Choose the appropriate genes.
      */

      i=0;

      g1Complete=FALSE;
      g2Complete=FALSE;

      i=0;

      while(! g1Complete)
      {
	 total+=augmentedParentScores[i];

	 if(g1 < total)
	 {
	    g1Complete=TRUE;
	    g1=i;
	    break;
	 }

	 i++;
      }

      total=0;
      i=0;

      while(! g2Complete)
      {
	 total+=augmentedParentScores[i];

	 if(g2 < total)
	 {
	    g2Complete=TRUE;
	    g2=i;
	    break;
	 }

	 i++;
      }

      /* Select a random position in the gene string */

      positionInGene = (int)(get_random_number_double(randomNumberSeed++) * lengthOfGene);

      /* Create the initial child gene */

      memcpy(childGenes[numberOfChildGenes],parentGenes[g1],positionInGene);
      childGenes[numberOfChildGenes][positionInGene]='\0';

      p=parentGenes[g2]+positionInGene;

      strcat(childGenes[numberOfChildGenes],p);

      /* Mutate the child gene by choosing a random position and mutating it at random */

      while(positionInGene >= lengthOfGene)
      {
	 positionInGene=(int)(get_random_number_double(randomNumberSeed++) * lengthOfGene);
      }

      randomNumber=get_random_number_double(randomNumberSeed++);

      if(randomNumber > PC)
      {
	 if(childGenes[numberOfChildGenes][positionInGene] == '1')
	 {
	    childGenes[numberOfChildGenes][positionInGene] = '0';
         }
         else
         {
            childGenes[numberOfChildGenes][positionInGene] = '1';
         }
      }

      /* Check if the number of 1's is greater than the threshold */

      if( number_of_occurrences(childGenes[numberOfChildGenes],'1') <= numberOfONPositions )
      {
	 numberOfChildGenes++;
	 continue;
      }

      /* Some positions have to be changed from 1 to 0 */

      while( number_of_occurrences(childGenes[numberOfChildGenes],'1') > numberOfONPositions )
      {
	 positionInGene=0;

	 while( childGenes[numberOfChildGenes][positionInGene] != '1' )
	 {
	    positionInGene=(int)(get_random_number_double(randomNumberSeed++) * lengthOfGene);
	 }

	 childGenes[numberOfChildGenes][positionInGene]='0';
      }

      numberOfChildGenes++;

   } /* End of while(numberOfChildGenes < population) */

} /* End of function "create_genes" */


double select_best_genes(char **parentGenes,double *parentScores,
		         char **childGenes,double *childScores,
		         int population,int lengthOfGene,
			 char *bestGene)
{
   double bestScore=0;

   int overallGeneIndex[2 * MAX_POPULATION];

   char bestGenes[2 * MAX_POPULATION][MAX_NUMBER_OF_INTERFACE_POSITIONS];

   double bestScores[2 * MAX_POPULATION];

   int i=0;

   /* Copy all the scores (child and parent scores) into a new array */

   for(i=0;i<population;i++)
   {
      overallScores[i]=parentScores[i];
      overallGeneIndex[i]=i;
   }

   for(i=0;i<population;i++)
   {
      overallScores[i+population] = childScores[i];
      overallGeneIndex[i+population] = i + population;
   }

   /* Arrange the scores in descending order using qsort */

   qsort((void *)overallGeneIndex, (size_t) 2*population, (size_t) sizeof(int), &compare);

   /* Now, choose the best ranking population genes */

   for(i=0;i<population;i++)
   {
      if(overallGeneIndex[i] >= population)
      {
	 strcpy(bestGenes[i],childGenes[overallGeneIndex[i]-population]);
	 bestScores[i]=childScores[overallGeneIndex[i]-population];
      }
      else
      {
	 strcpy(bestGenes[i],parentGenes[overallGeneIndex[i]]);
	 bestScores[i]=parentScores[overallGeneIndex[i]];
      }
   }

   /* Copy the best genes into the parent set */

   for(i=0;i<population;i++)
   {
      strcpy(parentGenes[i],bestGenes[i]);
      parentScores[i]=bestScores[i];
   }

   bestScore=overallScores[overallGeneIndex[0]];
   strcpy(bestGene,parentGenes[0]);

   /* Return the best score */

   return bestScore;

} /* End of function "select_best_genes" */


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-pdb") )
      {
	 strcpy(pdbListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-pos") )
      {
	 strcpy(interfacePositionsListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-ang") )
      {
	 strcpy(interfaceAnglesFilename,param[i+1]);
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

   pdbList=(char **)malloc(1000 * sizeof(char *));

   numberOfPDB=read_pdb_list(pdbListFilename,pdbList);

   /* Create the folds for training and validation of the network -> COMPLETED */

   folds=(char **)malloc(NUMBER_OF_FOLDS * sizeof(char *));

   create_fold_ranges(folds,numberOfPDB,NUMBER_OF_FOLDS);

   /* Read interface positions into a list -> Available in the library */

   interfacePositionsList=(char **)malloc(MAX_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char *));
   numberOfInterfacePositions=read_interface_positions_list(interfacePositionsList,interfacePositionsListFilename);

   /* Read the interface angles */

   interfaceAngles=(double *)malloc(numberOfPDB * sizeof(double));

   fp=fopen(interfaceAnglesFilename,"r");

   for(i=0;i<numberOfPDB;i++)
   {
      rewind(fp);

      interfaceAngles[i]=extract_interface_angle(fp,pdbList[i]);
   }

   fclose(fp);

   /* Create an array of linked lists for the Kabat numbering of all antibodies -> COMPLETED */

   kabatNumberingFirstPointers=(KABATLIST **)malloc(MAX_NUMBER_OF_PDB * sizeof(KABATLIST *));

   create_kabat_numbering_lists(pdbList,numberOfPDB,kabatNumberingFirstPointers);

   /* Allocate space for the parent and child genes */

   parentGenes=(char **)malloc(POPULATION * sizeof(char *));
   childGenes=(char **)malloc(POPULATION * sizeof(char *));

   for(i=0;i<POPULATION;i++)
   {
      parentGenes[i]=(char *)malloc((numberOfInterfacePositions + 1) * sizeof(char));
      childGenes[i]=(char *)malloc((numberOfInterfacePositions + 1) * sizeof(char));
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

   initialise_population(parentGenes,POPULATION,numberOfInterfacePositions,NUMBER_OF_ON_POSITIONS);

   /* While overall correlation coefficient of the gene set is less than 0.75, do the following steps */

   bestScore=0;

   parentScores=(double *)malloc(POPULATION * sizeof(double));
   childScores=(double *)malloc(POPULATION * sizeof(double));

   bestGene=(char *)malloc( (numberOfInterfacePositions + 1) * sizeof(char));

   while(bestScore < 0.75)
   {
      /* First, score the parent genes */

      i=0;

      while(i < POPULATION)
      {
	 parentScores[i]=calculate_coefficient_over_folds(parentGenes[i],
							  folds,
							  interfacePositionsList,
							  numberOfInterfacePositions,
							  &snnsFilenames);
	 i++;
      }

      /* Create a child gene set with POPULATION genes */

      create_genes(parentGenes,
		   parentScores,
		   childGenes,
		   POPULATION,
		   numberOfInterfacePositions,
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
