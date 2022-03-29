# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <strings.h>
# include <unistd.h>
# include <math.h>

# include "grid_ga.h"

static char inputFilename[MAX_FILENAME_LENGTH],
	    outputFilename[MAX_FILENAME_LENGTH];

static double thresholdScore = 0,
	      score = 0;

static int generation = 0;

static int numberOfGenesCommitted = 0,
	   geneNumber = 0;

static char *p = NULL,
	    line[2000],
	    ignore[100];

static char gene[MAX_NUMBER_OF_INTERFACE_POSITIONS + 1],
	    genes[10000][MAX_NUMBER_OF_INTERFACE_POSITIONS + 1];


/* ------------------------- FUNCTION DECLARATION SECTION --------------------- */

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

int read_genes_in_file(char *filename,char genes[][MAX_NUMBER_OF_INTERFACE_POSITIONS + 1]);

BOOL is_gene_redundant(char *gene,char genes[][MAX_NUMBER_OF_INTERFACE_POSITIONS + 1],int numberOfGenesCommitted);

/*---------------------- END OF FUNCTION DECLARATION SECTION ----------------- */


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>",argv[0]);
   printf("\n\nArguments are:\n");
   printf("\n1. -in <Input file>");
   printf("\n2. -out <Output file> - Optional");
   printf("\n3. -gen <Generation number to write into file> - Optional");
   printf("\n4. -th <Threshold score for selection of a gene>");
   printf("\n\n");

} /* End of function "Usage" */


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   BOOL inputFilenameFlag = FALSE,
	thresholdScoreFlag = FALSE;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-in") )
      {
	 inputFilenameFlag = TRUE;
	 strcpy(inputFilename,param[i+1]);
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
	 generation = atoi(param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-th") )
      {
	 thresholdScoreFlag = TRUE;
	 thresholdScore = atof(param[i+1]);
	 i+=2;
	 continue;
      }
      else
      {
	 return FALSE;
      }
   }

   if( (! thresholdScoreFlag) || (! inputFilenameFlag) )
   {
      return FALSE;
   }
   else
   {
      return TRUE;
   }

} /* End of function "parse_command_line_parameters" */


int read_genes_in_file(char *filename,char genes[][MAX_NUMBER_OF_INTERFACE_POSITIONS + 1])
{
   FILE *fp;

   char line[2000],
	ignore[100];

   int numberOfGenesCommitted = 0;

   fp = fopen(filename,"r");

   while( fgets(line,2000,fp) )
   {
      if( strstr(line,"GENE #") )
      {
	 sscanf(line,"%s%s%s",ignore,ignore,gene);

	 if(! is_gene_redundant(gene,genes,numberOfGenesCommitted) )
	 {
	    strcpy(genes[numberOfGenesCommitted],gene);
	    numberOfGenesCommitted++;
	 }
      }
   }

   fclose(fp);

   return numberOfGenesCommitted;

} /* End of function "read_genes_in_file" */


BOOL is_gene_redundant(char *gene,char genes[][MAX_NUMBER_OF_INTERFACE_POSITIONS + 1],int numberOfGenesCommitted)
{
   int i=0;

   while(i < numberOfGenesCommitted)
   {
      if(! strcmp(gene,genes[i]) )
      {
	 /* gene and genes[i] are identical. Return TRUE */

	 return TRUE;
      }

      i++;
   }

   return FALSE;

} /* End of function "is_gene_redundant" */


int main(int argc,char **argv)
{
   FILE *fp = NULL,
	*wfp = NULL;

   /* Parse command line parameters */

   outputFilename[0] = 0;

   if(argc < 4)
   {
      Usage(argv);
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   /* Check for file permissions and existence */

   if( access(inputFilename,R_OK) )
   {
      fprintf(stderr,"\nFile \"%s\" cannot be opened in read mode.\nAborting program\n\n",inputFilename);
      return 0;
   }

   /* Read input file and write output accordingly */

   if(outputFilename[0] == 0)
   {
      wfp = stdout;
      numberOfGenesCommitted = 0;
   }
   else
   {
      /* An argument for output file has been passed to the function. Check if the file
	 already exists and if it does, read the genes in the file.
      */

      if(! access(outputFilename,R_OK) )
      {
	 numberOfGenesCommitted = read_genes_in_file(outputFilename,genes);
      }
      else
      {
	 numberOfGenesCommitted = 0;
      }

      wfp = fopen(outputFilename,"a");
   }

   /* Read current input file and process genes in it */

   fp = fopen(inputFilename,"r");

   geneNumber = numberOfGenesCommitted + 1;

   while( fgets(line,5000,fp) )
   {
      /* Remove newline character at end of line */

      p = strchr(line,'\n');
      (*p) = '\0';

      if( strstr(line,"GENE #") )
      {
	 p=strchr(line,':'); /* GENE #1: 000010010101 */

	 p+=2;
	 strcpy(gene,p);

	 continue;
      }
      else
      if( strstr(line,"SCORE") )
      {
	 sscanf(line,"%s%lf",ignore,&score);

	 if(score >= thresholdScore)
	 {
	    /* Write gene and score into file.

	       Use the function "do_write_gene_into_file":

	       int do_write_gene_into_file(FILE *wfp,
	                                   int geneNumber,
	                                   char *gene,
	                                   double score,
	                                   int generation);
	    */

	    if( is_gene_redundant(gene,genes,numberOfGenesCommitted) )
	    {
	       continue;
	    }

	    do_write_gene_into_file(wfp,
				    geneNumber,
				    gene,
				    score,
				    generation);

	    strcpy(genes[geneNumber-1],gene);

	    geneNumber++;
	    numberOfGenesCommitted++;
	 }

	 continue;
      }
   }

   fclose(fp);
   fclose(wfp);

   /* Return 1 to the operating system */

   return 1;

} /* End of program */
