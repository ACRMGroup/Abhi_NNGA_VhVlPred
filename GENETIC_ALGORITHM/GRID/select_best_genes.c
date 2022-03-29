# define MAIN

# include <stdio.h>
# include <string.h>

# include "grid_ga.h"
# include "ga.h"

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

static int POPULATION=0,
	   i=0;

static char parentGenesFilename[200],
	    childGenesFilename[200];

static FILE *fp=NULL,
	    *wfp=NULL;

static struct gene_score_generation parentGeneScoreGen[MAX_POPULATION],
				    childGeneScoreGen[MAX_POPULATION];

static char bestGene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];

static char outputFilename[200];

static int generation=0;

/* --------------------------------- FUNCTION DECLARATION SECTION --------------------------------- */

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

/* ------------------------------ END OF FUNCTION DECLARATION SECTION ----------------------------- */


/* -------------------------------- FUNCTION DEFINITION SECTION ----------------------------------- */


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nArguments are:\n");
   printf("\n1. -mom <File with parent genes and scores>");
   printf("\n2. -kid <File with child genes and scores>");
   printf("\n3. -out <Output file to write best genes and their scores>");
   printf("\n4. -gen <Generation>");
   printf("\n\n");

} /* End of function "Usage" */

BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-mom") )
      {
	 strcpy(parentGenesFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-kid") )
      {
	 strcpy(childGenesFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(outputFilename,param[i+1]);
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
      {
	 return FALSE;
      }
   }

   return TRUE;

} /* End of function "parse_command_line_parameters" */

/* ----------------------------- END OF FUNCTION DEFINITION SECTION ------------------------------- */



int main(int argc,char **argv)
{
   if(argc < 8)
   {
      Usage(argv);
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      printf("\nError parsing command line parameters\n\n");
      Usage(argv);
      return 0;
   }

   /* Check if environment variable GENE_POPULATION has been set. Exit program if not */

   if(! getenv("GENE_POPULATION") )
   {
      printf("\nEnvironment variable \"GENE_POPULATION\" not set. Aborting program\n\n");
      return 0;
   }

   POPULATION=atoi(getenv("GENE_POPULATION"));

   if(POPULATION <= 0)
   {
      printf("\nEnvironment variable \"GENE_POPULATION\" not set properly. Aborting program\n\n");
      return 0;
   }

   /* Read parent genes with their scores. This is done using the function:

      void do_read_gene_from_file(FILE *fp,int geneNumber,struct gene_score_generation *geneScoreGen);

      This function reads the information associated with a gene into a structure of the form:

      struct gene_score_generation
      {
         char gene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];
         double score;
         int generation;
      };
   */

   GENE.parentGenes=(char **)malloc(POPULATION * sizeof(char *));
   GENE.parentScores=(double *)malloc(POPULATION * sizeof(double));

   fp=fopen(parentGenesFilename,"r");

   for(i=0;i<POPULATION;i++)
   {
      do_read_gene_from_file(fp,i+1,&(parentGeneScoreGen[i]));

      GENE.parentGenes[i] = (char *)malloc( (strlen(parentGeneScoreGen[i].gene) + 1) * sizeof(char));

      strcpy(GENE.parentGenes[i],parentGeneScoreGen[i].gene);

      GENE.parentScores[i]=parentGeneScoreGen[i].score;
   }

   fclose(fp);

   GENE.lengthOfGene=strlen(parentGeneScoreGen[0].gene);

   /* Read the child genes and their scores in a similar manner */

   GENE.childGenes=(char **)malloc(POPULATION * sizeof(char *));
   GENE.childScores=(double *)malloc(POPULATION * sizeof(double));

   fp=fopen(childGenesFilename,"r");

   for(i=0;i < POPULATION;i++)
   {
      do_read_gene_from_file(fp,i+1,&(childGeneScoreGen[i]));

      GENE.childGenes[i] = (char *)malloc( (strlen(childGeneScoreGen[i].gene) + 1) * sizeof(char));

      strcpy(GENE.childGenes[i],childGeneScoreGen[i].gene);

      GENE.childScores[i]=childGeneScoreGen[i].score;
   }

   fclose(fp);

   /* Now, find the best genes using the function "select_best_genes". The prototype of the function is:

      double select_best_genes(char **parentGenes,double *parentScores,
                               char **childGenes,double *childScores,
                               int population,int lengthOfGene,
                               char *bestGene)
   */

   select_best_genes_elitist(GENE.parentGenes,GENE.parentScores,
			     GENE.childGenes,GENE.childScores,
			     POPULATION,GENE.lengthOfGene,
			     bestGene);

   /* Write the best genes and their scores into the output file.

      Use the function "do_write_gene_into_file":

      int do_write_gene_into_file(FILE *wfp,
                                  int geneNumber,
                                  char *gene,
                                  double score,
                                  int generation);
   */

   wfp=fopen(outputFilename,"w");

   for(i=0;i<POPULATION;i++)
   {
      do_write_gene_into_file(wfp,
			      i+1,
			      GENE.parentGenes[i],
			      GENE.parentScores[i],
			      generation);
      /*
      fprintf(wfp,"GENE #%d: %s\nSCORE: %f\nGENERATION: %d\n",i+1,GENE.parentGenes[i],GENE.parentScores[i],generation);
      fputs("-------------------------------------\n",wfp);

      free(GENE.parentGenes[i]);
      free(GENE.childGenes[i]);
      */
   }

   fclose(wfp);

   free(GENE.parentScores);
   free(GENE.childScores);

   return 1;
}
