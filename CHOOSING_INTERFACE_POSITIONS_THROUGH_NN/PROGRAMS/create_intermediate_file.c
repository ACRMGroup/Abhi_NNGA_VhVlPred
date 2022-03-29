# include <string.h>
# include <stdio.h>
# include <stdlib.h>
# include <ctype.h>
# include "bioplib/SysDefs.h"

# define MAXIMUM_NUMBER_OF_INTERFACE_POSITIONS 20

/* ------------------------- GLOBAL VARIABLES DECLARATION SECTION ------------------------- */


typedef struct kabat_list
{
   char position[8],
	residueOneLetterCode;

   struct kabat_list *next;

}KABATLIST;

/* Variables associated with file names. */

static char interfacePositionsListFilename[100],
	    kabatNumberingFilename[100],
	    residuePropertiesMatrixFilename[100],
	    interfaceAngleFilename[100];

static char pdbCode[8];

static char outputFilename[100];

/* --------------------- END OF GLOBAL VARIABLES DECLARATION SECTION ---------------------- */


int read_interface_positions_list(char **interfacePositionsList,char *interfacePositionsListFilename)
{
   FILE *fp=NULL;
   char line[200];
   int i=0;

   fp=fopen(interfacePositionsListFilename,"r");

   if(! fp)
      return -1;

   strcpy(line,"");

   while( (line[0] != 'H') && (line[0] != 'L') )
      fgets(line,200,fp);

   do
   {
      if( (line[0] != 'L') && (line[0] != 'H') )
      {
	 fgets(line,200,fp);
	 continue;
      }

      line[strlen(line)-1]='\0';

      interfacePositionsList[i]=(char *)malloc(8 * sizeof(char *));

      sscanf(line,"%s",interfacePositionsList[i++]);

   } while( fgets(line,200,fp) );

   fclose(fp);

   return i; /* Return number of interface positions. */

} /* End of function "read_interface_positionslist". */


KABATLIST *read_kabat_numbered_file(char *kabatNumberingFilename)
{
   FILE *fp=NULL;

   KABATLIST *first=NULL,
	     *prev=NULL,
	     *p=NULL;

   char line[200];

   fp=fopen(kabatNumberingFilename,"r");

   if(! fp)
      return NULL;

   while( fgets(line,200,fp) )
   {
      p=(KABATLIST *)malloc(sizeof(KABATLIST));

      sscanf(line,"%s",p->position);

      p->residueOneLetterCode=line[strlen(line)-2];

      if(! first)
	 first=p;

      if(prev)
	 prev->next=p;

      prev=p;
   }

   fclose(fp);

   return first;

} /* End of function "read_kabat_numbered_file". */


BOOL find_residues_at_interface_positions(char **interfacePositionsList, 	/* Strings of interface positions. */
					  int numberOfInterfacePositions,	/* Number of elements in the array above. */
                                          KABATLIST *kabatNumberedListFirstPtr, /* Pointer to the first node of Kabat numbered list. */
					  char *residuesAtInterfacePositions)   /* Array to which residues must be written. */
{
   int i=0;
   KABATLIST *p=NULL;

   while(i < numberOfInterfacePositions)
   {
      p=kabatNumberedListFirstPtr;

      while( strcmp(p->position,interfacePositionsList[i]) && p )
	 p=p->next;

      if(! p)
	 return FALSE;

      residuesAtInterfacePositions[i]=p->residueOneLetterCode;

      i++;
   }

   return TRUE;

} /* End of function "find_residues_at_interface_positions". */


int correlate_interface_positions_and_matrix(char *residuePropertiesMatrixFilename, /* File with AA properties. */
                                             char *residuesAtInterfacePositions,    /* List of residues at interface positions. */
                                             int numberOfInterfacePositions,	    /* Number of elements in the array above. */
                                             float **interfaceResiduesMatrix)	    /* Matrix of Amino acid properties. */
{
   FILE *fp=fopen(residuePropertiesMatrixFilename,"r"),
	*start=NULL;

   int i=0,
       rewind=0;

   char line[200],
	startLine[200],
	ignore;

   if(! fp)
      return -1;

   strcpy(line,"");
   strcpy(startLine,"");

   while(! isalpha(line[0]) && fp )
      fgets(line,200,fp);

   strcpy(startLine,line);

   start=fp;

   while(i < numberOfInterfacePositions)
   {
      fseek(fp,-1*rewind,SEEK_CUR);

      rewind=0;
      strcpy(line,startLine);

      while(line[0] != residuesAtInterfacePositions[i])
      {
	 fgets(line,200,fp);
	 rewind+=strlen(line);
      }

      /* Once we are here, the required numbers are contained in the string "line". */

      interfaceResiduesMatrix[i]=(float *)malloc(4 * sizeof(float));

      sscanf(line,"%c %f %f %f %f",&ignore,
				   &interfaceResiduesMatrix[i][0],
				   &interfaceResiduesMatrix[i][1],
				   &interfaceResiduesMatrix[i][2],
				   &interfaceResiduesMatrix[i][3]);
      i++;
   }

   fclose(fp);

   return i;

} /* End of function "correlate_interface_positions_and_matrix". */


float extract_interface_angle(char *interfaceAngleFilename,char *pdbCode)
{
   FILE *fp=NULL;

   char line[200],
	ignore[30];

   float interfaceAngle=-1;

   fp=fopen(interfaceAngleFilename,"r");

   if(! fp)
      return -1;

   strcpy(line,"");

   while(! strstr(line,pdbCode) )
      fgets(line,200,fp);

   /* The line following the line containing the PDB code contains the interface angle */

   fgets(line,200,fp);

   sscanf(line,"%s %s %f",ignore,ignore,&interfaceAngle);

   return interfaceAngle;
}


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nThe arguments are:\n");
   printf("\n1. -ilist <File with list of Interface positions>");
   printf("\n2. -kabnum <File with the Kabat numbered residues>");
   printf("\n3. -mat <File with Matrix representation of the amino acid properties>");
   printf("\n4. -angle <File with interface angles for every antibody>");
   printf("\n5. -pdb <PDB Code of antibody>");
   printf("\n6. -out <Optional parameter - Output file>");
   printf("\n\n");

} /* End of function "Usage". */


BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   while(i < numberOfParam)
   {
      if(! strcmp(param[i],"-ilist") )
      {
	 strcpy(interfacePositionsListFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-kabnum") )
      {
	 strcpy(kabatNumberingFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-mat") )
      {
	 strcpy(residuePropertiesMatrixFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-angle") )
      {
	 strcpy(interfaceAngleFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-pdb") )
      {
	 strcpy(pdbCode,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcmp(param[i],"-out") )
      {
	 strcpy(outputFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
	 return FALSE;
   }

   return TRUE;

} /* End of function "parse_command_line_parameters". */


int main(int argc,char **argv)
{
   /* Step 1: Declare variables local to the main section of the code. */

   KABATLIST *kabatNumberedListFirstPtr=NULL,
	     *p=NULL,*prev=NULL; /* Pointers used to release the memory allocated to the linked list. */

   char **interfacePositionsList=NULL;

   int numberOfInterfacePositions=0,
       i=0;

   char *residuesAtInterfacePositions=NULL;

   float **interfaceResiduesMatrix=NULL;

   float interfaceAngle=-1;

   FILE *wfp=NULL;

   /* Step 2: Read all the command line parameters into program variables */

   if(argc < 11)
   {
      Usage(argv);
      exit(0);
   }

   strcpy(outputFilename,"");

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      exit(0);
   }

   /* Step 3: Allocate memory from heap area for the arrays. */

   interfacePositionsList=(char **)malloc(MAXIMUM_NUMBER_OF_INTERFACE_POSITIONS * sizeof(char **));

   /* Step 4: Read the interface positions into an array. */

   numberOfInterfacePositions=read_interface_positions_list(interfacePositionsList,interfacePositionsListFilename);

   /* Step 5: Read the Kabat numbered residues from the file into a linked list. */

   kabatNumberedListFirstPtr=read_kabat_numbered_file(kabatNumberingFilename);

   /* Step 6: Locate the residues corresponding to the interface positions. */

   residuesAtInterfacePositions=(char *)malloc(numberOfInterfacePositions * sizeof(char));

   find_residues_at_interface_positions(interfacePositionsList,
					numberOfInterfacePositions,
					kabatNumberedListFirstPtr,
					residuesAtInterfacePositions);

   /* Step 7: Make correspondence between the amino acids in the interface positions and
	      the matrix that represents their various properties.
   */

   interfaceResiduesMatrix=(float **)malloc(numberOfInterfacePositions * sizeof(float *));

   correlate_interface_positions_and_matrix(residuePropertiesMatrixFilename,
					    residuesAtInterfacePositions,
					    numberOfInterfacePositions,
					    interfaceResiduesMatrix);

   /* Step 8: Obtain the interface interface angle from the respective file */

   interfaceAngle=extract_interface_angle(interfaceAngleFilename,pdbCode);

   /* Step 9: Write the information gathered from the above steps into the output file. */

   if(! strcmp(outputFilename,"") )
   {
      wfp=stdout;
   }
   else
   {
      wfp=fopen(outputFilename,"a");
   }

   fprintf(wfp,"# PDB Code: %s\n\n",pdbCode);

   for(i=0;i<numberOfInterfacePositions;i++)
   {
      fprintf(wfp,"%s %c: %3.4f %3.4f %3.4f %3.4f\n",interfacePositionsList[i],
						     residuesAtInterfacePositions[i],
						     interfaceResiduesMatrix[i][0],
						     interfaceResiduesMatrix[i][1],
						     interfaceResiduesMatrix[i][2],
						     interfaceResiduesMatrix[i][3]);
   }

   fprintf(wfp,"\nInterface Angle: %3.4f\n----------------\n",interfaceAngle);

   fclose(wfp);

   /* Step 10: Release memory allocated to the variables */

   for(i=0;i<numberOfInterfacePositions;i++)
   {
      free(interfacePositionsList[i]);
      free(interfaceResiduesMatrix[i]);
   }

   free(interfacePositionsList);
   free(interfaceResiduesMatrix);

   p=kabatNumberedListFirstPtr;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   free(residuesAtInterfacePositions);

   return 0;
}
