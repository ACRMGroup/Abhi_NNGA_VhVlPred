# define MAIN 1

# include <stdio.h>
# include <string.h>
# include "ga.h"
# include "grid_ga.h"
# include "intermediate.h"

char outputFilename[200],
     parentGenesFilename[200],
     interfacePositionsListFilename[200];

int NUMBER_OF_ON_POSITIONS=20;

char **interfacePositionsList=NULL;

int numberOfInterfacePositions=0;

int i=0;

FILE *wfp=NULL;

int POPULATION=0;

/* -------------------------- FUNCTION DECLARATIONS ----------------------------- */

/* ga.h

   int compare(const void *a, const void *b);

   int number_of_occurrences(char *string,char ch);

   void initialise_population(struct gene *genes,       // Structure containing genes
                              int population,           // Number of genes
                              int lengthOfGene,         // Length of gene
                              int numberOfONPositions); // Number of positions to switch on (the character 1)

   int read_pdb_list(struct list_of_pdb *listOfPDB);

   void create_fold_ranges(char **folds,int numberOfPDB,int numberOfFolds);

   BOOL create_kabat_numbering_lists(struct list_of_pdb *listOfPDB,
                                     int numberOfPDB,
                                     KABATLIST **kabatNumberingFirstPointers);

   double calculate_coefficient_over_folds(char *gene,
                                           char **folds,
                                           int numberOfFolds,
                                           INTERFACE *interface,
                                           struct residue_properties *residueProperties,
                                           struct snns_filenames *snnsFilenames,
                                           KABATLIST **kabatNumberingFirstPointers,
                                           struct list_of_pdb *listOfPDB);

   double calculate_coefficient_over_single_fold(char *gene,
                                                 char *fold,
                                                 INTERFACE *interface,
                                                 struct residue_properties *residueProperties,
                                                 struct snns_filenames *snnsFilenames,
                                                 KABATLIST **kabatNumberingFirstPointers,
                                                 struct list_of_pdb *listOfPDB);

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
*/


void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

/* -------------------- END OF FUNCTION DECLARATION SECTION --------------------- */


/* ------------------------ FUNCTION DEFINITION SECTION ------------------------- */


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nArguments are:\n");
   printf("\n1. -int <File with list of interface positions>");
   printf("\n2. -out <Name of output file into which genes must be written (Default: parent_genes.txt)>");
   printf("\n3. -on <Maximum number of interface positions (Default: 20)>");
   printf("\n\n");
}

BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-int") )
      {
	 strcpy(interfacePositionsListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(parentGenesFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-on") )
      {
	 NUMBER_OF_ON_POSITIONS = atoi(param[i+1]);
	 i+=2;
	 continue;
      }
      else
      {
	 return FALSE;
      }
   }

   if(NUMBER_OF_ON_POSITIONS > 0)
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }

} /* End of function "parse_command_line_parameters" */


/* -------------------- END OF FUNCTION DEFINTION SECTION -------------------- */


int main(int argc,char **argv)
{
   /* Declare an object of type "struct gene".

      struct gene
      {
         char **parentGenes;
         char **childGenes;
         char *bestGene;

         int lengthOfGene;
         double *parentScores;
         double *childScores;
      };
   */

   struct gene genes;

   /* Parse command line parameters */

   parentGenesFilename[0] = 0;

   if(argc < 2)
   {
      Usage(argv);
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   if(! getenv("GENE_POPULATION") )
   {
      printf("\nPlease set environment variable \"GENE_POPULATION\" to the number of genes to be generated\n\n");
      return 0;
   }

   POPULATION=atoi(getenv("GENE_POPULATION"));

   if(parentGenesFilename[0] == 0)
   {
      strcpy(parentGenesFilename,"parent_genes.txt");
   }

   /* Read the list of interface positions from file. The function that does this is:

      int read_interface_positions_list(char **interfacePositionsList,char *interfacePositionsListFilename);
   */

   interfacePositionsList=(char **)malloc(MAX_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char *));

   numberOfInterfacePositions = read_interface_positions_list(interfacePositionsList,interfacePositionsListFilename);

   /* Now, create the parent genes. Function prototype:

      void initialise_population(struct gene *genes,       // Structure containing genes
                                 int population,           // Number of genes
                                 int numberOfONPositions); // Number of positions to switch on (the character 1)

      The strucure "gene" is defined in the following way:

      struct gene
      {
         char **parentGenes;
         char **childGenes;
         char *bestGene;
      
         int lengthOfGene;
         double *parentScores;
         double *childScores;
      };

      Before doing this, set the random number seed.
   */

   randomNumberSeed=getpid()+11;

   genes.lengthOfGene=numberOfInterfacePositions;

   genes.parentGenes=(char **)malloc(POPULATION * sizeof(char *));

   for(i=0;i<POPULATION;i++)
   {
      genes.parentGenes[i]=(char *)malloc( (genes.lengthOfGene + 1) * sizeof(char));
   }

   initialise_population(&genes,POPULATION,NUMBER_OF_ON_POSITIONS);

   /* Write the genes into the output file. The format of the output file is:

      GENE #1: Gene string

      SCORE: Whatever the score is, 0 if no score has been assigned as yet

      GENERATION: Generation - For parent genes, this is 1.

      Use the function "do_write_gene_into_file":

      int do_write_gene_into_file(FILE *wfp,
                                  int geneNumber,
                                  char *gene,
                                  double score,
                                  int generation)
   */

   wfp=fopen(parentGenesFilename,"w");

   for(i=0;i<POPULATION;i++)
   {
      do_write_gene_into_file(wfp,i+1,genes.parentGenes[i],0,1);

      /*

      fprintf(wfp,"GENE #%d: %s\nSCORE: 0\nGENERATION: 1\n",i+1,genes.parentGenes[i]);
      fputs("-------------------------------------\n",wfp);
      free(genes.parentGenes[i]);

      */
   }

   free(genes.parentGenes);

   fclose(wfp);

   return 1;

} /* End of program */
