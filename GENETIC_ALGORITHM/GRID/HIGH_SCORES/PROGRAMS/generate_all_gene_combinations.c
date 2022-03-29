# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include <strings.h>

# include "standard.h"
# include "bioplib/SysDefs.h"
# define MAX 100

static int numberOfBits=0,
	   generation=0,
	   numberOfFiles=0;

static char filenamePrefix[100];

static char outputFilename[MAX_FILENAME_LENGTH];

static FILE *wfp = NULL;

/* FUNCTION DECLARATION SECTION */

void reverseIt(char *temp,int j,FILE *wfp);

void printTheString(int i,int k,FILE *wfp);

void printGenes(int numberOfBits,int generation,int numberOfFiles);

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);

/* END OF FUNCTION DECLARATION SECTION */


void reverseIt(char *temp,int j,FILE *wfp)
{
   int i,k;
   char l;
   i = 0;
   k = j - 1;
  
   while(i < k)
   {
      l = temp[i];
      temp[i] = temp[k];
      temp[k] = l;
      i++;
      k--;            
   }

   temp[j] = '\0';
   fprintf(wfp,"%s\n",temp);
}

void printTheString(int i,int k,FILE *wfp)
{
   char str[MAX];
   int remainder;
   int j = 0; 
 
   while(i > 1)
   {
      remainder = i % 2;
      i = i/2;
 
      if(remainder == 0)
      {
	 str[j] = '0';
      }
      else
      {
	 str[j] = '1';
      }
    
      j++; 
   }

   
   if(i == 0) 
   {
      str[j] = '0';
   }
   else
   {
      str[j] = '1';
   }

   j++;

   while( j < k)
   { 
      str[j] = '0';
      j++;
   }

   str[j] = '\0';

   reverseIt(str,j,wfp); 
}

void printGenes(int numberOfBits,int generation,int numberOfFiles)
{
   long endNo = (long)(pow( (double)2, (double)numberOfBits) - 1);
   long i = 0;
   long numberPerFile = (long)(endNo/numberOfFiles);
   long outputFileNumber = 0;
   long int remainder=0;

   while(i <= endNo)
   {
      remainder = i % numberPerFile;

      if( remainder == 0 )
      {
	 outputFileNumber++;
	 sprintf(outputFilename,"%s_%ld.txt",filenamePrefix,outputFileNumber);

	 if(wfp)
	 {
	    fclose(wfp);
	 }

	 wfp = fopen(outputFilename,"w");
      }

      fprintf(wfp,"GENE #%ld: ",remainder + 1);

      printTheString(i,numberOfBits,wfp);

      fprintf(wfp,"SCORE: 0\nGENERATION: %d",generation);

      fputs("\n------------------------\n",wfp);

      i++; 
   }

   if(wfp)
   {
      fclose(wfp);
   }
}


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>",argv[0]);
   printf("\n\nArguments are:\n");
   printf("\n1. -bit <Number of bits>");
   printf("\n2. -gen <Generation to start with>");
   printf("\n3. -num <Number of files to be created> - Optional parameter (Default: 1)");
   printf("\n4. -pre <Prefix for file names>");
   printf("\n\n");
}


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-bit") )
      {
	 numberOfBits = atoi(param[i+1]);
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
      if(! strcasecmp(param[i],"-num") )
      {
	 numberOfFiles = atoi(param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-pre") )
      {
	 strcpy(filenamePrefix,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      {
	 return FALSE;
      }
   }

   if(numberOfFiles <= 0)
   {
      numberOfFiles = 1;
   }

   return TRUE;
}


int main(int argc,char **argv)
{
   if(argc < 3)
   {
      Usage(argv);
      return 0;
   }

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   if(numberOfBits > 32)
   {
      printf("\nThis program wont work for %d bits",numberOfBits);
      printf("\n\n");
      return 0;
   }

   printGenes(numberOfBits,generation,numberOfFiles);

   return 1;
}
