# include <stdio.h>
# include <string.h>
# include <strings.h>
# include <unistd.h>
# include <stdlib.h>

# include "gsl/gsl_statistics_double.h"
# include "arrays.h"
# include "bioplib/SysDefs.h"

static double scoresInGeneration[MAX_POPULATION];

static char inputFilename[MAX_FILENAME_LENGTH],
            outputFilename[MAX_FILENAME_LENGTH];

static char geneStrings[MAX_POPULATION][MAX_NUMBER_OF_INTERFACE_POSITIONS];

static int generation=0,
           numberOfRedundantSequences=0;

static double bestScoreInGeneration=0,
              averageScoreInGeneration=0,
              standardDeviationInGeneration=0,
              averageSimilarityInGeneration=0;

static BOOL printHeaderFlag = FALSE,
            printRedundancy = FALSE;

/* -------------------------- FUNCTION DEFINITION SECTION ------------------------ */

double degree_of_similarity(int g1,int g2)
{
   int numberOfCommonONPositions=0,
       totalNumberOfONPositions=0;

   int i=0,
       lengthOfGene =  strlen(geneStrings[g1]);

   while(i < lengthOfGene)
   {
      if( (geneStrings[g1][i] == '1') && (geneStrings[g2][i] == '1') )
      {
         numberOfCommonONPositions++;
         totalNumberOfONPositions+=2;

         i++;
         continue;
      }

      if(geneStrings[g1][i] == '1')
      {
         totalNumberOfONPositions++;
      }
      else
      if(geneStrings[g2][i] == '1')
      {
         totalNumberOfONPositions++;
      }

      i++;
   }

   return ( numberOfCommonONPositions/(double)(totalNumberOfONPositions) );

} /* End of function "degree_of_similarity" */


int calculate_statistics(double *bestScoreInGeneration,
                         double *averageScoreInGeneration,
                         double *standardDeviation,
                         int *numberOfRedundantSequences,
                         double *averageSimilarityInGeneration,
                         int numberOfGenesInGeneration)
{
   int i=0,
       j=0;

   int index[MAX_POPULATION];

   int numberOfUniqueSequences = 0;

   BOOL foundAMatchFlag = FALSE;

   double totalSimilarityInGene=0,
          totalSimilarityInGeneration=0;

   /* Calculate best score and average score */

   (*bestScoreInGeneration) = gsl_stats_max(scoresInGeneration, 1,numberOfGenesInGeneration);

   (*averageScoreInGeneration) = gsl_stats_mean(scoresInGeneration, 1, numberOfGenesInGeneration);

   (*standardDeviation) = gsl_stats_sd(scoresInGeneration, 1, numberOfGenesInGeneration);

   /* Find number of unique sequences */

   for(i=0;i < numberOfGenesInGeneration; i++)
   {
      /* Compare gene with all the members of the unique set */

      j=0;

      foundAMatchFlag = FALSE;

      while(j < numberOfUniqueSequences)
      {
         if(! strcmp(geneStrings[i],geneStrings[index[j]]) )
         {
            /* String geneStrings[i] is already present in the list of unique genes */

            foundAMatchFlag = TRUE;
            break;
         }

         j++;
      }

      if(foundAMatchFlag == FALSE)
      {
         /* geneStrings[i] doesnt have a presence in the list of unique genes.
            Add it to the list
         */

         index[numberOfUniqueSequences++] = i;
      }
   }

   (*numberOfRedundantSequences) = (numberOfGenesInGeneration - numberOfUniqueSequences);

   /* Find degree of redundancy for every sequence in the population */

   if(! printRedundancy)
   {
      return 1;
   }

   for(i=0;i < numberOfGenesInGeneration;i++)
   {
      totalSimilarityInGene=0;

      for(j=0;j < numberOfGenesInGeneration;j++)
      {
         totalSimilarityInGene+=degree_of_similarity(i,j);
      }

      totalSimilarityInGene = totalSimilarityInGene - 1; /* Similarity of the gene with itself */

      totalSimilarityInGeneration+=totalSimilarityInGene/(numberOfGenesInGeneration - 1);
   }

   (*averageSimilarityInGeneration) = totalSimilarityInGeneration/numberOfGenesInGeneration;

   /* Return from routine */

   return 1;

} /* End of function "calculate_statistics" */


int print_statistics(FILE *wfp,
                     int generation,
                     int numberOfRedundantSequences,
                     double bestScoreInGeneration,
                     double averageScoreInGeneration,
                     double standardDeviationInGeneration,
                     int numberOfGenesInGeneration,
                     double averageSimilarityInGeneration)
{
   fprintf(wfp,"%5d",generation);
   fprintf(wfp,"%15d",numberOfRedundantSequences);
   fprintf(wfp,"%19f",bestScoreInGeneration);
   fprintf(wfp,"%25f",averageScoreInGeneration);
   fprintf(wfp,"%23f",standardDeviationInGeneration);

   if(printRedundancy)
   {
      fprintf(wfp,"%30f",averageSimilarityInGeneration);
   }

   fprintf(wfp,"\n");

   return 1;
}


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>",argv[0]);
   printf("\n\nArguments are:\n");
   printf("\n1. -in <Input file>");
   printf("\n2. -out <Output file> - Optional");
   printf("\n3. -ph - Print header - Optional");
   printf("\n3. -r - Print redundancy - Optional");
   printf("\n\n");
}


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-in") )
      {
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
      if(! strcasecmp(param[i],"-ph") )
      {
         printHeaderFlag = TRUE;
         i++;
         continue;
      }
      if(! strcasecmp(param[i],"-r") )
      {
         printRedundancy = TRUE;
         i++;
         continue;
      }
      else
      {
         return FALSE;
      }
   }

   return TRUE;

} /* End of function "parse_command_line_parameters" */



int main(int argc,char **argv)
{
   FILE *fp,*wfp;

   char line[2000];
   char *p;

   int numberOfGenesInGeneration = -1,
       i = 0;

   fp = NULL;
   wfp = NULL;
   p = NULL;

   /* Parse command line parameters */

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

   /* Check for files */

   if( access(inputFilename,R_OK) )
   {
      printf("\nFile \"%s\" does not exist.\nAborting program\n\n",inputFilename);
      return 0;
   }

   if( access(outputFilename,W_OK) )
   {
      wfp = stdout;
   }
   else
   {
      wfp = fopen(outputFilename,"w");
   }

   /* Print the header for display */

   if(printHeaderFlag == TRUE)
   {
      fprintf(wfp,"\n");
      fprintf(wfp,"GENERATION");
      fprintf(wfp,"\tNO. OF");
      fprintf(wfp,"\t      BEST SCORE");
      fprintf(wfp,"\t      AVERAGE SCORE");
      fprintf(wfp,"\t      STANDARD");

      if(printRedundancy)
      {
         fprintf(wfp,"\t\t      AVERAGE SIMILARITY");
         fprintf(wfp,"\n\t       RED SEQS\t\t\t\t\t\t\t      DEVIATION\t\t          IN GENERATION\n");

         for(i=0;i<120;i++)
         {
            fputc('-',wfp);
         }

         fputc('\n',wfp);
      }
      else
      {
         fprintf(wfp,"\n\t       RED SEQS\t\t\t\t\t\t\t      DEVIATION\n");

         for(i=0;i<89;i++)
         {
            fputc('-',wfp);
         }

         fputc('\n',wfp);
      }

   }

   /* Read the contents of the file, calculate the statistics and print the values */

   fp=fopen(inputFilename,"r");

   while( fgets(line,2000,fp) )
   {
      p=strchr(line,'\n');
      (*p)='\0';

      if( strstr(line,"GENERATION #") )
      {
         /* GENERATION #1 */

         if(numberOfGenesInGeneration != -1)
         {
            /* Calculate the statistics using function "calculate_statistics" */

            numberOfGenesInGeneration++;

            calculate_statistics(&bestScoreInGeneration,
                                 &averageScoreInGeneration,
                                 &standardDeviationInGeneration,
                                 &numberOfRedundantSequences,
                                 &averageSimilarityInGeneration,
                                 numberOfGenesInGeneration);

            print_statistics(wfp,
                             generation,
                             numberOfRedundantSequences,
                             bestScoreInGeneration,
                             averageScoreInGeneration,
                             standardDeviationInGeneration,
                             numberOfGenesInGeneration,
                             averageSimilarityInGeneration);

            bestScoreInGeneration = 0;
            averageScoreInGeneration = 0;
            numberOfRedundantSequences = 0;
            numberOfGenesInGeneration = -1;
         }

         p=strchr(line,'#');
         p++;
         generation=atoi(p);
      }
      else
      if( strstr(line,"GENE ") )
      {
         p = strchr(line,':');
         p+=2;
         strcpy(geneStrings[++numberOfGenesInGeneration],p);
      }
      else
      if( strstr(line,"SCORE:") )
      {
         p=strchr(line,':');
         p+=2;
         scoresInGeneration[numberOfGenesInGeneration] = atof(p);
      }
   }

   /* Calculate and print statistics for the last generation */

   calculate_statistics(&bestScoreInGeneration,
                        &averageScoreInGeneration,
                        &standardDeviationInGeneration,
                        &numberOfRedundantSequences,
                        &averageSimilarityInGeneration,
                        numberOfGenesInGeneration);

   /* Print the statistics:

      int print_statistics(FILE *wfp,
                           int generation,
                           int numberOfRedundantSequences,
                           double bestScoreInGeneration,
                           double averageScoreInGeneration,
                           double standardDeviationInGeneration,
                           int numberOfGenesInGeneration,
                           double averageSimilarityInGeneration)
   */


   print_statistics(wfp,
                    generation,
                    numberOfRedundantSequences,
                    bestScoreInGeneration,
                    averageScoreInGeneration,
                    standardDeviationInGeneration,
                    numberOfGenesInGeneration,
                    averageSimilarityInGeneration);

   /* Close the file pointers */

   fclose(fp);
   fclose(wfp);

   /* Return 1 to the calling function */

   return 1;
}
