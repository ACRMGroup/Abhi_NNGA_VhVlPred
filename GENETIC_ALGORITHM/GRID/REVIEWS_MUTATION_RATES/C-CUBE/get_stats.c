# include <stdio.h>
# include <string.h>
# include <strings.h>
# include <unistd.h>
# include <errno.h>
# include <stdlib.h>

# include "gsl/gsl_statistics_double.h"
# include "arrays.h"
# include "general.h"
# include "bioplib/SysDefs.h"
# include "bioplib/array.h"

# define MAX_LENGTH_OF_GENE_SEGMENT 24
# define MAX_NUMBER_OF_PARTS 20
# define MAXBUFFER 2000

static double scoresInGeneration[MAX_POPULATION];

static char inputFilename[MAX_FILENAME_LENGTH],
            outputFilename[MAX_FILENAME_LENGTH],
            stringToSearch[20];

long int **numberInGeneSegment = NULL;

static int generation=0,
           numberOfRedundantSequences=0,
           numberOfDivisionsInGene = 0;

static double bestScoreInGeneration=0,
              averageScoreInGeneration=0,
              standardDeviationInGeneration=0,
              averageSimilarityInGeneration=0;

static BOOL printHeaderFlag = FALSE,
            printRedundancy = FALSE;

static int startGeneration = -1;


/* -------------------------- FUNCTION DEFINITION SECTION ------------------------ */

int split_gene_into_numbers(char *gene, long int *numberInGeneSegment)
{
   int j = 0,
       terminate = 0,
       currentLB = -1,
       currentUB = -1;

   char geneSegment[MAX_LENGTH_OF_GENE_SEGMENT];

   while(! terminate)
   {
      currentLB = currentUB + 1;
      currentUB = currentLB + MAX_LENGTH_OF_GENE_SEGMENT - 1;

      if( currentUB >= strlen(gene) )
      {
         currentUB = strlen(gene) - 1;
         terminate = 1;
      }

      get_substring(gene, geneSegment, currentLB, currentUB);

      numberInGeneSegment[j++] = strtol(geneSegment, NULL, 2);
   }

   return j;
}

int count_number_of_bits(int number)
{
   int c = 0;

   for(; number ; c++)
   {
      number &= number - 1;
   }

   return c;

} /* End of function "count_number_of_bits" */


double degree_of_similarity(int g1, int g2, int numberOfDivisionsInGene)
{
   int numberOfBits1 = 0,
       numberOfBits2 = 0,
       numberOfCommonBits = 0,
       i = 0;

   double degreeOfSimilarity = 0;

   for(i = 0 ; i < numberOfDivisionsInGene ; i++)
   {
      numberOfBits1 += count_number_of_bits(numberInGeneSegment[g1][i]);
      numberOfBits2 += count_number_of_bits(numberInGeneSegment[g2][i]);
      numberOfCommonBits += count_number_of_bits(numberInGeneSegment[g1][i] & numberInGeneSegment[g2][i]);
   }

   degreeOfSimilarity = (double)numberOfCommonBits/(numberOfBits1 + numberOfBits2);

   return degreeOfSimilarity;

} /* End of function "degree_of_similarity". */



int calculate_statistics(double *bestScoreInGeneration,
                         double *averageScoreInGeneration,
                         double *standardDeviation,
                         int *numberOfRedundantSequences,
                         double *averageSimilarityInGeneration,
                         int numberOfGenesInGeneration)
{
   int i=0,
       j=0,
       k=0,
       xorValue=0;

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

   for(i=0 ; i < numberOfGenesInGeneration ; i++)
   {
      /* Compare gene with all the members of the unique set */

      j = 0;

      foundAMatchFlag = FALSE;

      for(j = 0 ; j < numberOfUniqueSequences ; j++)
      {
         k = 0;
         xorValue = 0;

         while(k < numberOfDivisionsInGene)
         {
            xorValue += numberInGeneSegment[i][k] ^ numberInGeneSegment[index[j]][k];
            k++;
         }

         if(xorValue == 0)
         {
            foundAMatchFlag = TRUE;
            break;
         }
         else
         {
            foundAMatchFlag = FALSE;
         }
      }

      if(foundAMatchFlag == FALSE)
      {
         /* numberInGeneSegment[i] doesn't have a presence in the list of unique genes.
            Add it to the list.
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

   for(i = 0 ; i < numberOfGenesInGeneration ; i++)
   {
      totalSimilarityInGene=0;

      for(j=0;j < numberOfGenesInGeneration;j++)
      {
         totalSimilarityInGene+=degree_of_similarity(i, j, numberOfDivisionsInGene);
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
   printf("\n4. -gen <Generation> - Perform analysis on a certain generation - Optional");
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
      else
      if(! strcasecmp(param[i],"-r") )
      {
         printRedundancy = TRUE;
         i++;
         continue;
      }
      else
      if(! strcasecmp(param[i],"-gen") )
      {
         startGeneration = atoi(param[i + 1]);
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



int main(int argc,char **argv)
{
   FILE *fp,*wfp;

   char line[MAXBUFFER];
   char *p = NULL,
        eraseOption = 0;

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

   if(! access(outputFilename, W_OK) )
   {
      fprintf(stderr, "\nFile \"%s\" already exists. Erase? (Y/N)", outputFilename);
      scanf("%c", &eraseOption);

      if( (eraseOption == 'Y') || (eraseOption == 'y') )
      {
         wfp = fopen(outputFilename, "w");
      }
      else
      if( (eraseOption == 'N') || (eraseOption == 'n') )
      {
         wfp = stdout;
      }
      else
      {
         return 0;
      }
   }
   else
   {
      wfp = stdout;
   }

   /* Allocate memory for the genes to hold integer values of the genes */

   numberInGeneSegment = (long int **)Array2D(sizeof(long int), MAX_POPULATION, MAX_NUMBER_OF_PARTS);

   /* Check if a certain generation number is required. */

   fp=fopen(inputFilename,"r");

   if(startGeneration > 0)
   {
      sprintf(stringToSearch, " #%d: ", startGeneration);

      strcpy(line, "");

      while( (! strstr(line, stringToSearch) ) && ( ! feof(fp) ) )
      {
         fgets(line, MAXBUFFER, fp);
      }

      if( feof(fp) )
      {
         fclose(fp);

         fprintf(stderr, "\nGeneration #%d not found in file.\nAborting program.\n\n", startGeneration);

         FreeArray2D((char **)numberInGeneSegment, MAX_POPULATION, MAX_NUMBER_OF_PARTS);

         return(0);
      }

      /* Rewind the file pointer by the length of the string */

      fseek(fp, -1 * strlen(line), SEEK_CUR);
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

   numberOfGenesInGeneration = 0;

   while( fgets(line,MAXBUFFER,fp) )
   {
      p=strchr(line,'\n');
      (*p)='\0';

      if( strstr(line,"GENERATION #") )
      {
         /* GENERATION #1 */

         if(numberOfGenesInGeneration > 0)
         {
            /* Calculate the statistics using function "calculate_statistics" */

            /* numberOfGenesInGeneration++; */

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
            numberOfGenesInGeneration = 0;
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

         numberOfDivisionsInGene = split_gene_into_numbers(p, numberInGeneSegment[numberOfGenesInGeneration]);
         numberOfGenesInGeneration++;
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

   if(numberOfGenesInGeneration > 0)
   {
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
   }

   /* Release all the memory */

   FreeArray2D((char **)numberInGeneSegment, MAX_POPULATION, MAX_NUMBER_OF_PARTS);

   /* Close the file pointers */

   fclose(fp);
   fclose(wfp);

   /* Return 1 to the calling function */

   return 1;
}
