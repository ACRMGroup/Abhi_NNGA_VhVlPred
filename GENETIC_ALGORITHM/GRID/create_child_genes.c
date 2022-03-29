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

/* Format of struct gene_score_generation:

   struct gene_score_generation
   {
      char gene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];
      double score;
      int generation;
   };
*/

static struct gene_score_generation geneScoreGen[MAX_POPULATION];

/* Format of struct gene:

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

static struct gene GENE;

static int POPULATION=0;

static FILE *fp=NULL,
	    *wfp=NULL;

static char parentGenesFilename[200],
	    childGenesFilename[200];

static int numberOfONPositions=20;

static int i=0,
	   generation=0;

static char selectionProcedure[100];

/* ---------------------------- FUNCTION DECLARATION SECTION -------------------------- */

void free_memory();

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

/* ------------------------ END OF FUNCTION DECLARATION SECTION ----------------------- */


/* ---------------------------- FUNCTION DEFINITION SECTION --------------------------- */

void free_memory()
{
   int i=0;

   /* GENE.parentGenes - number of elements to be freed is given by POPULATION */

   for(i=0;i<POPULATION;i++)
   {
      free(GENE.parentGenes[i]);
   }

   free(GENE.parentGenes);

   /* GENE.parentScores */

   free(GENE.parentScores);

   /* GENE.childGenes - number of elements to be freed is given by POPULATION */

   for(i=0;i<POPULATION;i++)
   {
      free(GENE.childGenes[i]);
   }

} /* End of function "free_memory" */


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nArguments are:\n");
   printf("\n1. -in <Name of file containing parent genes and their associated scores>");
   printf("\n2. -gen <Generation>");
   printf("\n3. -out <Output file for child genes (Default: child_genes.txt)");
   printf("\n4. -proc <Selection procedure for the selection of parent genes>");
   printf("\n\n");

} /* End of function "Usage" */


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-in") )
      {
	 strcpy(parentGenesFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(childGenesFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-gen") )
      {
	 generation=atoi(param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-proc") )
      {
	 strcpy(selectionProcedure,param[i+1]);
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


/* ------------------------ END OF FUNCTION DEFINITION SECTION ------------------------ */


int main(int argc,char **argv)
{
   if(argc < 6)
   {
      Usage(argv);
      return 0;
   }

   selectionProcedure[0]=0;

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   /* Check if selection procedure has been indicated. If not, report error and exit program */

   if(selectionProcedure[0] == 0)
   {
      printf("\nSelection procedure not input. Aborting program\n\n");
      return 0;
   }

   /* Check if the environment variable GENE_POPULATION has been set. If not, report error and exit */

   if(! getenv("GENE_POPULATION") )
   {
      printf("\nEnvironment variable \"GENE_POPULATION\" not set. Exiting program\n\n");
      return 0;
   }
   else
   {
      POPULATION = atoi(getenv("GENE_POPULATION"));
   }

   if(POPULATION <= 0)
   {
      printf("\nEnvironment variable \"GENE_POPULATION\" not set properly. Exiting program\n\n");
      return 0;
   }

   /* Read the genes in the input file and scores. This is done using the function:

      void do_read_gene_from_file(FILE *fp,int geneNumber,struct gene_score_generation *geneScoreGen);
   */

   fp=fopen(parentGenesFilename,"r");

   if(! fp)
   {
      printf("\nFile \"%s\" does not exist. Exiting program\n\n",parentGenesFilename);
      return 0;
   }

   GENE.parentGenes=(char **)malloc(POPULATION * sizeof(char *));
   GENE.parentScores=(double *)malloc(POPULATION * sizeof(double));

   for(i=0;i < POPULATION;i++)
   {
      do_read_gene_from_file(fp,i+1,&(geneScoreGen[i]));

      GENE.parentGenes[i]=(char *)malloc( (strlen(geneScoreGen[i].gene) + 1) * sizeof(char));
      GENE.parentScores[i]=geneScoreGen[i].score;

      strcpy(GENE.parentGenes[i],geneScoreGen[i].gene);
   }

   fclose(fp);

   GENE.lengthOfGene=strlen(geneScoreGen[0].gene);

   /* Now, create the child genes. The function that does this is:

      void create_child_genes(struct gene *GENE,
			      int population,
			      int numberOfONPositions);
   */

   GENE.childGenes=(char **)malloc(POPULATION * sizeof(char *));

   for(i=0;i<POPULATION;i++)
   {
      GENE.childGenes[i]=(char *)malloc( (GENE.lengthOfGene + 1) * sizeof(char));
   }

   create_child_genes(&GENE,
		      POPULATION,
		      numberOfONPositions,
		      selectionProcedure);

   if(childGenesFilename[0] == 0)
   {
      strcpy(childGenesFilename,"child_genes.txt");
   }

   wfp=fopen(childGenesFilename,"w");

   /* Write all the child genes into a file.

      Use the function "do_write_gene_into_file":

      int do_write_gene_into_file(FILE *wfp,
                                  int geneNumber,
                                  char *gene,
                                  double score,
                                  int generation);
   */

   for(i=0;i<POPULATION;i++)
   {
      do_write_gene_into_file(wfp,
			      i+1,
			      GENE.childGenes[i],
			      0,
			      generation);
      /*

      fprintf(wfp,"GENE #%d: %s\nSCORE: 0\nGENERATION: %d\n",i+1,GENE.childGenes[i],generation);
      fputs("-------------------------------------\n",wfp);

      */
   }

   fclose(wfp);

   /* Release memory allocated to the variables and exit from program */

   free_memory();

   return 1;

} /* End of program */
