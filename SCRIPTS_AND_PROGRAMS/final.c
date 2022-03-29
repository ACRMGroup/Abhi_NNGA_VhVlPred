# include <stdio.h>
# include <stdarg.h>
# include <stdlib.h>
# include <string.h>
# include <strings.h>
# include <time.h>
# include <ctype.h>
# include <unistd.h>
# include <math.h>

# include "bioplib/SysDefs.h"

# define CORRELATION_PROGRAM_PATH "/home/bsm/martin/bin/correlation"
# define TEMP_PATH "/tmp/SNNS"
# define NUMBER_OF_SNNS_TRAINING_PARAMETERS 7

/* --------------------------- INTRODUCTION ---------------------------------


 This is an all-step encompassing program which goes through the following steps:

 1. Default values:
		    -> init function - Randomize weights - Randomize_Weights (parameter to setInitFunc).
		    -> learning function - Resilient propogation - Rprop (parameter to setLearnFunc).
		    -> update function - Topological order - Topological_Order (parameter to setUpdateFunc).
		    -> pruning function - Magnitude based pruning - MagPruning (parameter to setPruningFunc).
		    -> shuffling mode - TRUE - ( setShuffle(TRUE) ).

 2. Accept weights for training network. If nothing is input, use default weights.

			-> setInitFunc - setInitFunc(Randomize_Weights, 1.0, -1.0)

			-> setLearnFunc - setLearnFunc(Rprop,0.2,0,1e-6).

			-> setUpdateFunc - setUpdateFunc(Topological_Order).

			-> setPruningFunc - setPruningFunc("MagPruning", "Rprop", 15.0, 3.5, FALSE, 500, 90,1e6, 1.0)

			-> setShuffle - setShuffle(TRUE).

 3. Accept type of untrained network:
					-> 20 input layers, 10 hidden layers.
					-> 20 input layers, 20 hidden layers.
					-> 20 input layers, 30 hidden layers.

 4. Inputs for training set boundary - the remaining patterns automatically become validation patterns.

 5. Train the neural network and save the network.

 6. Validate with the validation set.

 7. Write out predicted and actual results into a file.

 8. Calculate pearson coefficient and error rate for the network.


------------------------- END OF INTRODUCTION ------------------------------ */


/* --------------- GLOBAL VARIABLES DECLARATION SECTION -------------------- */

static struct snns_filenames
{
   /* Parameters pertaining to filenames */

   char snnsTrainingScriptFilename[1000],
	snnsValidationScriptFilename[1000],
	untrainedNetworkFilename[1000],
	trainingPatternsFilename[1000],
	validationPatternsFilename[1000],
	trainedNetworkFilename[1000],
	predictedValuesOutputFilename[1000],
	inputOutputComparisonsFilename[1000];

}snnsFilenames;

static struct snns_parameters
{
   /* Parameters for training the network */

   int numberOfCycles;

   float sseThreshold;

   /* Functions and associated parameters */

   /* setInitFunc(Randomize_Weights, 1.0, -1.0) */

   char initFunction[20];
   float initFunctionParam1,initFunctionParam2;

   /* setLearnFunc(Rprop,0.2,0,1e-6) */

   char learningFunction[20];
   float learningFunctionParam1,
	 learningFunctionParam2;
   char learningFunctionParam3[8];

   /* setUpdateFunc(Topological_Order) */

   char updateFunction[20];

   /* setPruningFunc("MagPruning", "Rprop", 15.0, 3.5, FALSE, 500, 90,1e6, 1.0) */

   char pruningFunction[20];
   char pruningFunctionParam1[20];
   float pruningFunctionParam2,pruningFunctionParam3;
   char pruningFunctionParam4[8];
   float pruningFunctionParam5,pruningFunctionParam6;
   char pruningFunctionParam7[8];
   float pruningFunctionParam8;

   /* setShuffle(TRUE) */

   char shuffling[8];

   /* Priority flag for over riding number of cycles or vice versa or using a combination
      of user specified inputs for both
   */

   int ssePriority;
   int cyclesPriority;

}snnsParameters;


typedef struct snns_pattern_entry
{
   char pdbCode[8];
   char aminoAcidParameters[200][8];
   float interfaceAngle;

   struct snns_pattern_entry *next;

}node;


static int NUMBER_OF_RESIDUE_PARAMETERS=4;

static int userNumberOfCycles = -1;

static float userSSEThreshold = -1;

static char snnsParametersConfigFilename[100],
	    coefficientOutputFilename[100];

static char *abhiDataDir=NULL;

static int numberOfInputUnits=0,
           numberOfOutputUnits=1;

static char intermediateFilename[100];

static char trainingSetBoundaryString[20];

static int numberOfTrainingPatterns=0,
           numberOfValidationPatterns=0;

static node *trainingStartPointer=NULL,
            *validationStartPointer=NULL;

static int currentProcessID=-1;

FILE *wfp=NULL;

static float maximumInterfaceAngle=-10000,
	     minimumInterfaceAngle=1000;

static int numberOfPDB=0;


/* ----------- END OF GLOBAL VARIABLE DECLARATION SECTION ------------------ */


/* --------------- SUB - ROUTINES DECLARATION SECTION ---------------------- */

/*

int read_intermediate_file(char *intermediateFilename,
			   char *trainingSetBoundaryString,
			   node *trainingStartPointer,int *numberOfTrainingPatterns,
			   node *validationStartPointer,int *numberOfValidationPatterns);
*/

int read_intermediate_file(char *intermediateFilename,
			   char *trainingSetBoundaryString,
			   int *numberOfTrainingPatterns,
			   int *numberOfValidationPatterns);

BOOL split_boundary_string(char *trainingSetBoundaryString,
			   int *trainingBoundaryStart,int *trainingBoundaryEnd);

BOOL is_training_pattern(int currentPDBNumber,
                         int trainingStartPointer,int trainingEnd,
                         int numberOfPDB);

void write_patterns_into_file(FILE *wfp,
                              node *startPointer,
                              int numberOfPatterns,
                              int numberOfInputUnits,
                              int numberOfOutputUnits);

void write_snns_pattern_file_header(FILE *wfp,
                                    int numberOfPatterns,
                                    int numberOfInputUnits,
                                    int numberOfOutputUnits);

void write_snns_commands(struct snns_filenames *snnsFilenames,
                         struct snns_parameters *snnsParameters);

BOOL execute_script_and_calculate_correlation_coefficient(node *validationFirst,struct snns_filenames *snnsFilenames);

BOOL initialise_snns_parameters(char *configFilename,struct snns_parameters *snnsParameters);

void assign_snns_filenames(struct snns_filenames *snnsFilenames,int currentProcessID);

void free_list();

void remove_temporary_files();

void Usage(char **argv);

BOOL parse_command_line_parameters(int numberOfParam,char **param);


/* ------------ END OF SUB - ROUTINES DECLARATION SECTION ------------------ */



/* ---------------- SUB - ROUTINES DEFINITION SECTION ---------------------- */



/* int read_intermediate_file(char *intermediateFilename,
                              char *trainingSetBoundaryString,
                              node *trainingStartPointer,int *numberOfTrainingPatterns,
                              node *validationStartPointer,int *numberOfValidationPatterns):

   This function reads the intermediate file and based on the training set
   boundaries specified by the user, the appropriate patterns are copied
   into two linked lists of type "snns_pattern_entry".

   The function returns the number of input units that go into the neural
   network (which is in turn, the total number of residue parameters for
   every pattern).
*/

/*
int read_intermediate_file(char *intermediateFilename,
			   char *trainingSetBoundaryString,
			   node *trainingStartPointer,int *numberOfTrainingPatterns,
			   node *validationStartPointer,int *numberOfValidationPatterns)
*/
int read_intermediate_file(char *intermediateFilename,
			   char *trainingSetBoundaryString,
			   int *numberOfTrainingPatterns,
			   int *numberOfValidationPatterns)
{
   /* Step : Declare the variables to be used in the function */

   char line[500],ignore[8];

   char delimiters[]=" ",
	*token,*running; /* Variables for parsing the line into tokens */

   node *p=NULL,
	*prev=NULL,
	*first=NULL;

   node *trainingPrev=NULL,
	*validationPrev=NULL;

   int trainingBoundaryStart=-1,
       trainingBoundaryEnd=-1;

   int i=0,j=0;

   /* Step : Open the file and check if the pointer is NULL (if file does not exist).
	     In such a case, return -1 to the calling routine.
   */

   FILE *fp=fopen(intermediateFilename,"r");

   if(! fp)
   {
      return -1;
   }

   /* Step : Parse the training set boundary string. If it cant be parsed successfully,
	     return a value of -2 to the calling function.
   */

   if(! split_boundary_string(trainingSetBoundaryString,&trainingBoundaryStart,&trainingBoundaryEnd) )
   {
      return -2;
   }

   /* Step : Go through the intermediate file and isolate the following fields:

	     1. PDB Codes.
	     2. Amino acid parameters for every interface position.
	     3. Interface angles.

	     Create a linked list of these fields.
   */

   i=0;

   while( fgets(line,500,fp) )
   {
      line[strlen(line)-1] = '\0';

      if( strstr(line,"PDB") )
      {
	 /* Example of Line:

	    # PDB Code: 12e8
	 */

	 p=(node *)malloc(sizeof(node));

	 sscanf(line,"%s %s %s %s",ignore,ignore,ignore,p->pdbCode);

	 numberOfPDB++;

	 continue;
      }

      if( (line[0] == 'H') || (line[0] == 'L') )
      {
	 /* The number of parameters in the line is given by "NUMBER_OF_RESIDUE_PARAMETERS".
	    Example of Line:

	    H45 L: 4.0000 3.0000 0.5300 0.0000

	    The first two tokens are not essential, while the next set of parameters need to
	    be stored in an array. We use the function "strsep" to tokenize the string.
	 */

	 running=line;
	 token=NULL;

	 for(j=0;j < (NUMBER_OF_RESIDUE_PARAMETERS + 2);j++)
	 {
	    token=strsep(&running,delimiters);

	    if(j > 1)
	    {
	       strcpy(p->aminoAcidParameters[i++],token);
	    }
	 }

	 continue;
      }

      if( strstr(line,"Interface") )
      {
	 sscanf(line,"%s %s %f",ignore,ignore,&p->interfaceAngle);

	 if(p->interfaceAngle < 0)
	 {
	    p->interfaceAngle = -1 * p->interfaceAngle;
	 }

	 numberOfInputUnits=i;
	 i=0;

	 /* Scan for the greatest and the least interface angles */

	 if(minimumInterfaceAngle > p->interfaceAngle)
	 {
	    minimumInterfaceAngle = p->interfaceAngle;
	 }

	 if(maximumInterfaceAngle < p->interfaceAngle)
	 {
	    maximumInterfaceAngle = p->interfaceAngle;
	 }
      }

      if(! first)
      {
	 first=p;
	 p->next=NULL;
	 prev=p;
      }

      if(prev)
      {
	 prev->next=p;
	 prev=p;
      }
   }

   /* Step : Check if the boundaries for the training set are within the number of PDB.
	      If they arent, then return -3 to the calling routine.
   */

   if( (trainingBoundaryStart > numberOfPDB) || (trainingBoundaryEnd > numberOfPDB) )
   {
      return -3;
   }

   /* Step : Now that the creation of the linked list is complete, split the overall
	      list into 2 parts - training and validation sets based on the training
	      set boundaries.
   */

   p=first;
   i=1;

   *numberOfTrainingPatterns=0;
   *numberOfValidationPatterns=0;

   while(i <= numberOfPDB)
   {
      if( is_training_pattern(i,trainingBoundaryStart,trainingBoundaryEnd,numberOfPDB) )
      {
	 /* Move node to the training pattern list */

	 if(! trainingStartPointer)
	    trainingStartPointer=p;

	 if(trainingPrev)
	    trainingPrev->next=p;

	 trainingPrev=p;

	 (*numberOfTrainingPatterns)++;
      }
      else
      {
	 /* Move node to the validation pattern list */

	 if(! validationStartPointer)
	    validationStartPointer=p;

	 if(validationPrev)
	    validationPrev->next=p;

	 validationPrev=p;

	 (*numberOfValidationPatterns)++;
      }

      p=p->next;

      i++;

   } /* End of while loop. */

   /* Terminate the 2 lists */

   trainingPrev->next=NULL;
   validationPrev->next=NULL;

   /* Step : Return a value of 1 to the calling routine, indicating successful processing! */

   return numberOfInputUnits;

} /* End of function "read_intermediate_file" */


BOOL split_boundary_string(char *trainingSetBoundaryString,
                           int *trainingBoundaryStart,int *trainingBoundaryEnd)
{
   int i=0;

   if(! strchr(trainingSetBoundaryString,'-') )
   {
      return FALSE;
   }

   i=0;

   while( i < strlen(trainingSetBoundaryString) )
   {
      if( isdigit(trainingSetBoundaryString[i]) )
      {
	 i++;
	 continue;
      }

      if(trainingSetBoundaryString[i] == '-')
      {
	 i++;
	 continue;
      }

      return FALSE;
   }

   *trainingBoundaryStart=atoi(trainingSetBoundaryString);
   *trainingBoundaryEnd=atoi(strchr(trainingSetBoundaryString,'-')+1);

   return TRUE;

} /* End of function "split_boundary_string" */


BOOL is_training_pattern(int currentPDBNumber,
                         int trainingStartPointer,int trainingEnd,
                         int numberOfPDB)
{
   if(trainingStartPointer < trainingEnd)
   {
      /* If boundaries are: 1-400 */

      if( (currentPDBNumber >= trainingStartPointer) && (currentPDBNumber <= trainingEnd) )
      {
	 /* If current PDB number is 100 */

	 return TRUE;
      }
      else
      {
	 /* If current PDB number is 430 */

	 return FALSE;
      }
   }
   else
   {
      /* This feature is a little more complicated than the first case.

	 If there are 500 structures in the set and boundaries are: 400-200.

	 trainingStartPointer=400.
	 trainingEnd=200.

	 The boundaries may be rewritten are: 1-200 and 400-500.
      */

      if( (currentPDBNumber <= trainingEnd) ||
	  ( (currentPDBNumber >= trainingStartPointer) && (currentPDBNumber <= numberOfPDB) ) )
      {
	 /* If current PDB number is 100 */

	 return TRUE;
      }
      else
      {
	 /* If current PDB number is 430 */

	 return FALSE;
      }
   }

   return FALSE;

} /* End of function "is_training_pattern". */


/* void write_patterns_into_file(FILE *wfp,
                                 node *startPointer,
                                 int numberOfPatterns,
                                 int numberOfInputUnits,
                                 int numberOfOutputUnits)

   This program writes patterns from a linked list of type "node" into a file. A FILE pointer
   to the file in write mode is passed to the function as argument.
*/

void write_patterns_into_file(FILE *wfp,
			      node *startPointer,
			      int numberOfPatterns,
			      int numberOfInputUnits,
			      int numberOfOutputUnits)
{
   node *p;

   int i=0;

   float numerator=0,
	 denominator=(maximumInterfaceAngle - minimumInterfaceAngle);

   write_snns_pattern_file_header(wfp,numberOfPatterns,numberOfInputUnits,numberOfOutputUnits);

   p=startPointer;

   while(p)
   {
      fprintf(wfp,"# PDB Code: %s\n",p->pdbCode);

      for(i=0;i<numberOfInputUnits;i++)
      {
	 fprintf(wfp,"%s ",p->aminoAcidParameters[i]);

	 if( ( (i+1) % NUMBER_OF_RESIDUE_PARAMETERS) == 0 )
	 {
	    fputc('\n',wfp);
	 }
      }

      /* Before printing the interface angle, scale it to a value between 0 and 1.
	 This is done in the following way:

	 interface angle fraction (scaled) =    (interface angle - least interface angle)
					     ------------------------------------------------
					     (maximum interface angle - least interface angle)
      */

      numerator=(p->interfaceAngle - minimumInterfaceAngle);

      /* fprintf(wfp,"# Interface angle for PDB %s\n%f\n",p->pdbCode,p->interfaceAngle); */

       fprintf(wfp,"# Interface angle for PDB %s\n%f\n",p->pdbCode,(numerator/denominator));

      p=p->next;
   }

} /* End of function "write_patterns_into_file". */


/* void write_snns_pattern_file_header(FILE *wfp,
                                       int numberOfPatterns,
                                       int numberOfInputUnits,
                                       int numberOfOutputUnits):

   This function writes the first few lines of a pattern file, as required by SNNS.
*/

void write_snns_pattern_file_header(FILE *wfp,
				    int numberOfPatterns,
				    int numberOfInputUnits,
				    int numberOfOutputUnits)
{
   time_t t;
   char *str;

   time(&t);
   str = ctime(&t);

   fputs("SNNS pattern definition file V3.2\n",wfp);
   fprintf(wfp,"generated at %s\n\n\n",str);

   fprintf(wfp,"No. of patterns : %d\n",numberOfPatterns);
   fprintf(wfp,"No. of input units : %d\n",numberOfInputUnits);
   fprintf(wfp,"No. of output units : %d\n\n",numberOfOutputUnits);

} /* End of function "write_snns_pattern_file_header". */


/* void write_snns_commands(struct snns_filenames *snnsFilenames,
                            struct snns_parameters *snnsParameters)

   This function performs the task of writing the SNNS commands so that
   they may be executed using batchman. The function accepts a FILE
   pointer and a pointer to the structure snns_filenames as arguments.
*/

void write_snns_commands(struct snns_filenames *snnsFilenames,
			 struct snns_parameters *snnsParameters)
{
   /* Step 1: Declare variables local to the function */

   FILE *wfp=NULL;

   /* Step 2: Write commands to train the network. The values to be written into the
  	      commands file are stored in 2 structures:

      static struct snns_filenames
      {
         char snnsTrainingScriptFilename[1000],
	      snnsValidationScriptFilename[1000],
              untrainedNetworkFilename[1000],
              trainingPatternsFilename[1000],
              validationPatternsFilename[1000],
              trainedNetworkFilename[1000],
              predictedValuesOutputFilename[1000];

      }snnsFilenames;
      
      static struct snns_parameters
      {
         int numberOfCycles;
         
         float sseThreshold;
         
         // setInitFunc(Randomize_Weights, 1.0, -1.0)
         
         char initFunction[20];
         float initFunctionParam1,initFunctionParam2;
         
         // setLearnFunc - setLearnFunc(Rprop,0.2,0,1e-6).
         
         char learningFunction[20];
         float learningFunctionParam1,
               learningFunctionParam2,
         char learningFunctionParam3[8];
         
         // setUpdateFunc - setUpdateFunc(Topological_Order).
         
         char updateFunction[20];
         
         // setPruningFunc - setPruningFunc("MagPruning", "Rprop", 15.0, 3.5, FALSE, 500, 90,1e6, 1.0)
         
         char pruningFunction[20];
         char pruningFunctionParam1[20];
         float pruningFunctionParam2,pruningFunctionParam3;
         char pruningFunctionParam4[8];
         float pruningFunctionParam5,pruningFunctionParam6;
         char pruningFunctionParam7[8];
         float pruningFunctionParam8;
         
         // setShuffle - setShuffle(TRUE).
         
         char shuffling[8];

	 // Priority flag for over riding number of cycles using SSE threshold.

	 int ssePriority;
	 int cyclesPriority;
      
      }snnsParameters;
   */

   wfp=fopen(snnsFilenames->snnsTrainingScriptFilename,"w");

   /* loadNet("encoder.net") */

   fprintf(wfp,"loadNet(\"%s\")",snnsFilenames->untrainedNetworkFilename);
   fputs("\n\n",wfp);

   /* loadPattern("encoder.pat") */

   fprintf(wfp,"loadPattern(\"%s\")",snnsFilenames->trainingPatternsFilename);
   fputs("\n\n",wfp);

   /* setInitFunc("Randomize_Weights", 1.0, -1.0):

      char initFunction[20];
      float initFunctionParam1,initFunctionParam2;
   */

   fprintf(wfp,"setInitFunc(\"%s\",%f,%f)",snnsParameters->initFunction,
					   snnsParameters->initFunctionParam1,
					   snnsParameters->initFunctionParam2);
   fputs("\n\n",wfp);

   /* setLearnFunc(Rprop,0.2,0,1e-6):

      char learningFunction[20];
      float learningFunctionParam1,
            learningFunctionParam2,
      char learningFunctionParam3[8];
   */

   fprintf(wfp,"setLearnFunc(\"%s\",%f,%f,%s)",snnsParameters->learningFunction,
					       snnsParameters->learningFunctionParam1,
					       snnsParameters->learningFunctionParam2,
					       snnsParameters->learningFunctionParam3);
   fputs("\n\n",wfp);

   /* setUpdateFunc(Topological_Order):

      char updateFunction[20];
   */

   fprintf(wfp,"setUpdateFunc(\"%s\")",snnsParameters->updateFunction);
   fputs("\n\n",wfp);

   /* setPruningFunc("MagPruning", "Rprop", 15.0, 3.5, FALSE, 500, 90,1e6, 1.0)

      char pruningFunction[20];
      char pruningFunctionParam1[20];
      float pruningFunctionParam2,pruningFunctionParam3;
      char pruningFunctionParam4[8];
      float pruningFunctionParam5,pruningFunctionParam6;
      char pruningFunctionParam7[8];
      float pruningFunctionParam8;
   */

   fprintf(wfp,"setPruningFunc(\"%s\",\"%s\",%f,%f,%s,%f,%f,%s,%f)",snnsParameters->pruningFunction,
							            snnsParameters->pruningFunctionParam1,
								    snnsParameters->pruningFunctionParam2,
								    snnsParameters->pruningFunctionParam3,
								    snnsParameters->pruningFunctionParam4,
								    snnsParameters->pruningFunctionParam5,
								    snnsParameters->pruningFunctionParam6,
								    snnsParameters->pruningFunctionParam7,
								    snnsParameters->pruningFunctionParam8);
   fputs("\n\n",wfp);

   /* setShuffle(TRUE) */

   fprintf(wfp,"setShuffle(%s)",snnsParameters->shuffling);
   fputs("\n\n",wfp);

   /* initNet() */

   fputs("initNet()\n\n",wfp);

   /* while SSE > 6.9 and CYCLES < 1000 do
         trainNet()
      endwhile
   */

   if(snnsParameters->ssePriority && ! snnsParameters->cyclesPriority)
   {
      fprintf(wfp,"while SSE > %f do\n\n",snnsParameters->sseThreshold);
   }
   else
   if(snnsParameters->cyclesPriority && ! snnsParameters->ssePriority)
   {
      fprintf(wfp,"while CYCLES < %d do\n\n",snnsParameters->numberOfCycles);
   }
   else
   {
      fprintf(wfp,"while SSE > %f and CYCLES < %d do\n\n",snnsParameters->sseThreshold,
						          snnsParameters->numberOfCycles);
   }

   fputs("   trainNet()\n\nendwhile\n\n",wfp);

   /* saveNet("encoder.trained.net") */

   fprintf(wfp,"saveNet(\"%s\")\n\n",snnsFilenames->trainedNetworkFilename);

   /* print ("Cycles trained: ", CYCLES) */

   fputs("print(\"Cycles trained: \",CYCLES)\n\n",wfp);

   /* print ("Training stopped at error: ", SSE) */

   fputs("print(\"Training stopped at error: \",SSE)\n\n",wfp);

   /* This concludes writing the training script. Close file pointer */

   fclose(wfp);

   /* Step 2: Write the script to validate the trained network */

   wfp=fopen(snnsFilenames->snnsValidationScriptFilename,"w");

   /* loadNet("trained_network.net") */

   fprintf(wfp,"loadNet(\"%s\")\n\n",snnsFilenames->trainedNetworkFilename);

   /* loadPattern("validation_pattern_file.pat") */

   fprintf(wfp,"loadPattern(\"%s\")\n\n",snnsFilenames->validationPatternsFilename);

   /* testNet() */

   fputs("testNet()\n\n",wfp);

   /* saveResult("validation_result.res") */

   fprintf(wfp,"saveResult(\"%s\")",snnsFilenames->predictedValuesOutputFilename);

   fclose(wfp);

} /* End of function "write_snns_commands" */


/* BOOL execute_script_and_calculate_correlation_coefficient(node *validationFirst,struct snns_filenames *snnsFilenames):

   This function correlates the actual interface angles and the interface angles
   predicted by the neural network. A seperate file that contains the actual and
   predicted interface angles is written out. A pearson coefficient calculation
   is performed against this file to correlate the accuracy of the interface angle
   predictions by the neural network.
*/

BOOL execute_script_and_calculate_correlation_coefficient(node *validationFirst,struct snns_filenames *snnsFilenames)
{
   /* Step 1: Declare the local variables required */

   FILE *fp=NULL,
	*wfp=NULL,
	*out=NULL;

   node *p=validationFirst;

   float actualInterfaceAngles[1000],
	 predictedInterfaceAngles[1000];

   float squaredError=0;

   char line[500],previousLine[500],command[200];

   int actualInterfaceAngleCount=0,
       predictedInterfaceAngleCount=0,
       i=0;

   FILE *pipe=NULL;

   float denominator=maximumInterfaceAngle - minimumInterfaceAngle;

   /* Step 2: Run batchman on the SNNS commands file */

   command[0]=0;

   sprintf(command,"batchman -q < %s;batchman -q < %s",snnsFilenames->snnsTrainingScriptFilename,
						       snnsFilenames->snnsValidationScriptFilename);
   /* system(command); */

   pipe=popen(command,"r");
   pclose(pipe);

   /* Step 3: Read the actual interface angles from the validation list */

   while(p)
   {
      actualInterfaceAngles[actualInterfaceAngleCount++] = (p->interfaceAngle - minimumInterfaceAngle)/denominator;
      p=p->next;
   }

   /* Step 4: Read the output file of the neural network to get the predicted angles */

   fp=fopen(snnsFilenames->predictedValuesOutputFilename,"r");

   line[0]=0;

   while( fgets(line,500,fp) && (line[0] != '#') ); /* Till line is: #1.1 */

   fgets(line,500,fp); /* Move to next line. Line is now the first set of input patterns: 4 3 0.53 0 10 6 0.37 0 3 2 */

   while( fgets(line,500,fp) )
   {
      line[strlen(line)-1] = '\0';

      if(line[0] == '#')
      {
	 /* Line before this contains the interface angle */

	 predictedInterfaceAngles[predictedInterfaceAngleCount++]=atof(previousLine);
      }

      strcpy(previousLine,line);
   }

   /* Add the last predicted interface angle to the array of predicted angles */

   predictedInterfaceAngles[predictedInterfaceAngleCount++]=atof(line);

   fclose(fp);

   /* Step 4: Compare the number of predicted and actual interface angles. If they dont match
	      return a FALSE to the calling routine
   */

   if(predictedInterfaceAngleCount != actualInterfaceAngleCount)
   {
      return FALSE;
   }

   /* Step 5: Now that the numbers are clear, write the two arrays into a file */

   if(coefficientOutputFilename[0] != 0)
   {
      out=fopen(coefficientOutputFilename,"a");
   }
   else
   {
      out=stdout;
   }

   fprintf(out,"Training set range: %s\n",trainingSetBoundaryString);

   if( (numberOfPDB - numberOfTrainingPatterns) > 1 )
   {
      wfp=fopen(snnsFilenames->inputOutputComparisonsFilename,"w");

      i=0;

      while(i < predictedInterfaceAngleCount)
      {
	 fprintf(wfp,"%f %f\n",actualInterfaceAngles[i],predictedInterfaceAngles[i]);

	 i++;
      }

      fclose(wfp);

      /* Calculate the pearsons correlation coefficient on the actual and predicted
	 interface angles using the program "correlation".
      */

      sprintf(command,"%s %s",CORRELATION_PROGRAM_PATH,snnsFilenames->inputOutputComparisonsFilename);

      pipe=popen(command,"r");

      fgets(line,200,pipe);

      fputs(line,out);

      pclose(pipe);
   }
   else
   {
      /* This is a case of boot strapping where all but one pattern have been used to train the network.
	 In this case, we calculate the error between the actual and predicted values and output it.
      */

      squaredError=pow((predictedInterfaceAngles[0] - actualInterfaceAngles[0]),2);

      fprintf(out,"Squared error and positive error: %4.4f %4.4f\n",squaredError,sqrt(squaredError));
   }

   fputs("################################################\n",out);

   fclose(out);

   return TRUE;

} /* End of function "execute_script_and_calculate_correlation_coefficient" */


/* BOOL initialise_snns_parameters(char *configFilename,struct snns_parameters *snnsParameters):

   This function initialises training parameters for the neural network. A file that contains
   the parameters along with suitable headers is received as an argument by the function. A
   TRUE is returned from the function following successful reading and parsing of the configuration
   file.
*/

BOOL initialise_snns_parameters(char *configFilename,struct snns_parameters *snnsParameters)
{
   /* Step 1: Declare all variables local to the function */

   FILE *fp=fopen(configFilename,"r");

   int fieldCount=0;

   char line[200],
	ignore[100];

   /* Step 2: Check if the configuration file exists. If not, return FALSE to the calling routine. */

   if(! fp)
   {
      return FALSE;
   }

   /* Step 3: Go through the remaining lines in the configuration file and record the values into the
	      appropriate variables in the structure.
   */

   while( fgets(line,200,fp) )
   {
      line[strlen(line)-1] = '\0';

      if( ! isupper(line[0]) )
      {
	 continue;
      }

      if( strstr(line,"CYCLES") )
      {
	 /* CYCLES: 200

	    int numberOfCycles;
	 */

	 sscanf(line,"%s %d",ignore,
			     &(snnsParameters->numberOfCycles));
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"SSE") )
      {
	 /* SSE: 1.5

	    float sseThreshold;
	 */

	 sscanf(line,"%s %f",ignore,
			     &(snnsParameters->sseThreshold));
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"INIT_FUNC:") )
      {
	 /* INIT_FUNC: Randomize_Weights 1.0 -1.0

            char initFunction[20];
            float initFunctionParam1,initFunctionParam2;
	 */

	 sscanf(line,"%s %s %f %f",ignore,
				   snnsParameters->initFunction,
				   &(snnsParameters->initFunctionParam1),
				   &(snnsParameters->initFunctionParam2));
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"LEARN_FUNC") )
      {
	 /* LEARN_FUNC: Rprop 0.1 50.0 1e-6

            char learningFunction[20];
            float learningFunctionParam1,
                  learningFunctionParam2,
            char learningFunctionParam3[8];
	 */

	 sscanf(line,"%s %s %f %f %s",ignore,
				      snnsParameters->learningFunction,
				      &(snnsParameters->learningFunctionParam1),
				      &(snnsParameters->learningFunctionParam2),
				      snnsParameters->learningFunctionParam3);
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"UPDATE_FUNC") )
      {
	 /* UPDATE_FUNC: Topological_Order

	    char updateFunction[20];
	 */

	 sscanf(line,"%s %s",ignore,
			     snnsParameters->updateFunction);
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"PRUNING_FUNC") )
      {
	 /* PRUNING_FUNC: MagPruning Rprop 15.0 3.5 FALSE 500 90 1e6 1.0

            char pruningFunction[20];
            char pruningFunctionParam1[20];
            float pruningFunctionParam2,pruningFunctionParam3;
            char pruningFunctionParam4[8];
            float pruningFunctionParam5,pruningFunctionParam6;
            char pruningFunctionParam7[8];
            float pruningFunctionParam8;
	 */

	 sscanf(line,"%s %s %s %f %f %s %f %f %s %f",ignore,
						     snnsParameters->pruningFunction,
						     snnsParameters->pruningFunctionParam1,
						     &(snnsParameters->pruningFunctionParam2),
						     &(snnsParameters->pruningFunctionParam3),
						     snnsParameters->pruningFunctionParam4,
						     &(snnsParameters->pruningFunctionParam5),
						     &(snnsParameters->pruningFunctionParam6),
						     snnsParameters->pruningFunctionParam7,
						     &(snnsParameters->pruningFunctionParam8));
	 fieldCount++;
	 continue;
      }

      if( strstr(line,"SHUFFLE_FLAG") )
      {
	 /* SHUFFLE_FLAG: TRUE

	    char shuffling[8];
	 */

	 sscanf(line,"%s %s",ignore,
			     snnsParameters->shuffling);
	 fieldCount++;
	 continue;
      }
   }

   snnsParameters->ssePriority=0;
   snnsParameters->cyclesPriority=0;

   if(fieldCount == NUMBER_OF_SNNS_TRAINING_PARAMETERS)
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }

} /* End of function "initialise_snns_parameters". */


/* void assign_snns_filenames(snns_filenames *snnsFilenames,int currentProcessID):

   This function is used to assign names to various files that are to be used in
   the program. The definition of the structure snns_filenames is as follows:

   static struct snns_filenames
   {
      char snnsTrainingScriptFilename[1000],
	   snnsValidationScriptFilename[1000],
           untrainedNetworkFilename[1000],
           trainingPatternsFilename[1000],
           validationPatternsFilename[1000],
           trainedNetworkFilename[1000],
           predictedValuesOutputFilename[1000],
           inputOutputComparisonsFilename[1000];
   };
*/

void assign_snns_filenames(struct snns_filenames *snnsFilenames,int currentProcessID)
{
   char *currentWorkingDirectory;

   currentWorkingDirectory=(char *)malloc(500 * sizeof(char));

   getcwd(currentWorkingDirectory,500);

   /* snnsTrainingScriptFilename: File that stores the commands for the SNNS batchman script */

   sprintf(snnsFilenames->snnsTrainingScriptFilename,"%s/%d_training.cmd",TEMP_PATH,currentProcessID);   

   /* snnsValidationScriptFilename: File that stores the commands for the SNNS batchman script */

   sprintf(snnsFilenames->snnsValidationScriptFilename,"%s/%d_validation.cmd",TEMP_PATH,currentProcessID);   

   /* untrainedNetworkFilename: The untrained Network path. This is supplied by the user */

   /* trainingPatternFilename: File into which training patterns must be written */

   sprintf(snnsFilenames->trainingPatternsFilename,"%s/%d_training.pat",TEMP_PATH,currentProcessID);

   if(snnsFilenames->trainingPatternsFilename[0] != '/')
   {
      sprintf(snnsFilenames->trainingPatternsFilename,"%s/%s",currentWorkingDirectory,
							      snnsFilenames->trainingPatternsFilename);
   }

   /* validationPatternFilename: File into which validation patterns must be written */

   sprintf(snnsFilenames->validationPatternsFilename,"%s/%d_validation.pat",TEMP_PATH,currentProcessID);

   if(snnsFilenames->validationPatternsFilename[0] != '/')
   {
      sprintf(snnsFilenames->validationPatternsFilename,"%s/%s",currentWorkingDirectory,
								snnsFilenames->validationPatternsFilename);
   }

   /* trainedNetworkFilename: File where trained network is saved */

   sprintf(snnsFilenames->trainedNetworkFilename,"%s/%d_trained.net",TEMP_PATH,currentProcessID);

   if(snnsFilenames->trainedNetworkFilename[0] != '/')
   {
      sprintf(snnsFilenames->trainedNetworkFilename,"%s/%s",currentWorkingDirectory,
							    snnsFilenames->trainedNetworkFilename);
   }

   /* predictedValuesOutputFilename: File that stores the neural network predictions for the validation patterns */

   sprintf(snnsFilenames->predictedValuesOutputFilename,"%s/%d.out",TEMP_PATH,currentProcessID);

   if(snnsFilenames->predictedValuesOutputFilename[0] != '/')
   {
      sprintf(snnsFilenames->predictedValuesOutputFilename,"%s/%s",currentWorkingDirectory,
								   snnsFilenames->predictedValuesOutputFilename);
   }

   /* inputOutputComparisonsFilename: File that stores the actual interface angle fractions along with
      those that have been predicted by the neural network.
   */

   sprintf(snnsFilenames->inputOutputComparisonsFilename,"%s/%d_io_comparison.txt",TEMP_PATH,currentProcessID);

   if(snnsFilenames->inputOutputComparisonsFilename[0] != '/')
   {
      sprintf(snnsFilenames->inputOutputComparisonsFilename,"%s/%s",currentWorkingDirectory,
								    snnsFilenames->inputOutputComparisonsFilename);
   }

   free(currentWorkingDirectory);

} /* End of function "assign_snns_filenames". */


/* void free_list(): This function releases memory used by the linked list of patterns */

void free_list()
{
   node *p=NULL,
	*prev=NULL;

   p=trainingStartPointer;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

   p=validationStartPointer;

   while(p)
   {
      prev=p;
      p=p->next;
      free(prev);
   }

} /* End of function "free_list". */


/* void remove_temporary_files(): This function removes all the temporary files created in the program */

void remove_temporary_files()
{
   /* In the structure snns_filenames, the following files have to be removed.

      static struct snns_filenames
      {
         char snnsTrainingScriptFilename[1000],		-> To be removed
              snnsValidationScriptFilename[1000],	-> To be removed
              untrainedNetworkFilename[1000],
              trainingPatternsFilename[1000],		-> To be removed
              validationPatternsFilename[1000],		-> To be removed
              trainedNetworkFilename[1000],		-> To be removed
              predictedValuesOutputFilename[1000],	-> To be removed
              inputOutputComparisonsFilename[1000];	-> To be removed.
      };
   */

   remove(snnsFilenames.snnsTrainingScriptFilename);
   remove(snnsFilenames.snnsValidationScriptFilename);
   remove(snnsFilenames.trainingPatternsFilename);
   remove(snnsFilenames.validationPatternsFilename);
   remove(snnsFilenames.trainedNetworkFilename);
   remove(snnsFilenames.predictedValuesOutputFilename);
   remove(snnsFilenames.inputOutputComparisonsFilename);

} /* End of function "remove_temporary_files" */


/* void Usage():

   This function lists the command line options of the program.
*/


void Usage(char **argv)
{
   printf("\nUsage: %s <Arguments>\n",argv[0]);
   printf("\nArguments are:\n");
   printf("\n1. -net <Path of untrained network file>");
   printf("\n3. -trb <Training set boundary string (Eg. 1-400)>");
   printf("\n3. -int <Path of intermediate file>");
   printf("\n4. -sse <Sum-of-Squares error threshold while training the network> - Optional parameter");
   printf("\n5. -scfg <File containing the SNNS configuration parameters> - Optional parameter");
   printf("\n6. -cyc <Number of cycles that the network must be trained on> - Optional parameter");
   printf("\n\n");

} /* End of function "Usage". */



BOOL parse_command_line_parameters(int numberOfParam,char **param)
{
   int i=1;

   int parameterCount=0,
       minimumNumberOfParameters=3;

   char *currentWorkingDirectory=NULL;

   currentWorkingDirectory=(char *)malloc(500 * sizeof(char));

   getcwd(currentWorkingDirectory,500);

   while(i < numberOfParam)
   {
      if(! strcasecmp(param[i],"-int") )
      {
	 strcpy(intermediateFilename,param[i+1]);
	 i+=2;
	 parameterCount++;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-scfg") )
      {
	 strcpy(snnsParametersConfigFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-net") )
      {
	 if(param[i+1][0] != '/')
	 {
	    /* The path is relative to the current working directory. Prefix the current working directory
	       to filename.
	    */

	    sprintf(snnsFilenames.untrainedNetworkFilename,"%s/%s",currentWorkingDirectory,param[i+1]);
	 }
	 else
	 {
	    strcpy(snnsFilenames.untrainedNetworkFilename,param[i+1]);
	 }

	 i+=2;
	 parameterCount++;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-sse") )
      {
	 userSSEThreshold=atof(param[i+1]);
	 snnsParameters.ssePriority=1;
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-cyc") )
      {
	 userNumberOfCycles=atoi(param[i+1]);
	 snnsParameters.cyclesPriority=1;
	 i+=2;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-trb") )
      {
	 strcpy(trainingSetBoundaryString,param[i+1]);
	 i+=2;
	 parameterCount++;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(coefficientOutputFilename,param[i+1]);
	 i+=2;
	 continue;
      }
      else
      {
	 return FALSE;
      }
   }


   free(currentWorkingDirectory);

   if(parameterCount != minimumNumberOfParameters)
   {
      return FALSE;
   }
   else
   {
      return TRUE;
   }

} /* End of function "parse_command_line_parameters". */


/* -------------------- MAIN SECTION OF THE CODE --------------------------- */

int main(int argc,char **argv)
{
   /* Step : Check if required number of inputs have been supplied by user. If not
	      print usage scheme of the program and exit.
   */

   if(argc < 6)
   {
      Usage(argv);
      return 0;
   }

   /* Step : Parse the command line parameters. If any problems of mismatch are
	      encountered, print usage scheme of program and exit.

	      Before entering the routine, make NULL, all parameters that are
	      optional in the program.
   */

   snnsParametersConfigFilename[0]=0;
   coefficientOutputFilename[0]=0;
   userSSEThreshold=-1;
   userNumberOfCycles=-1;

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   /*
   numberOfInputUnits=read_intermediate_file(intermediateFilename,
					     trainingSetBoundaryString,
					     trainingStartPointer,&numberOfTrainingPatterns,
					     validationStartPointer,&numberOfValidationPatterns);
   */

   numberOfInputUnits=read_intermediate_file(intermediateFilename,
					     trainingSetBoundaryString,
					     &numberOfTrainingPatterns,
					     &numberOfValidationPatterns);

   if(numberOfInputUnits < 0)
   {
      printf("\n\nError reading file \"%s\" or parsing boundary string \"%s\"\n\n",intermediateFilename,
									           trainingSetBoundaryString);
      return 0;
   }

   currentProcessID=getpid();

   /* Step : Set parameters for the SNNS file names and SNNS training parameters */

   assign_snns_filenames(&snnsFilenames,currentProcessID);

   if(snnsParametersConfigFilename[0] == 0)
   {
      abhiDataDir=getenv("ABHIDATADIR");
      sprintf(snnsParametersConfigFilename,"%s/ABNUM/snns_training_parameters.txt",abhiDataDir);
   }

   if(! initialise_snns_parameters(snnsParametersConfigFilename,&snnsParameters) )
   {
      printf("\nError reading file \"%s\" or file does not exist.\nAborting program\n\n",snnsParametersConfigFilename);
      return 0;
   }

   /* Over-ride default SSE threshold or number of training cycles if they have been specified by the user */

   if(userSSEThreshold != -1)
   {
      snnsParameters.sseThreshold=userSSEThreshold;
   }

   if(userNumberOfCycles != -1)
   {
      snnsParameters.numberOfCycles=userNumberOfCycles;
   }

   /* Step : Write the training and validation set patterns into respective files */

   wfp=fopen(snnsFilenames.trainingPatternsFilename,"w");

   write_patterns_into_file(wfp,
                            trainingStartPointer,
                            numberOfTrainingPatterns,
                            numberOfInputUnits,
                            numberOfOutputUnits);
   fclose(wfp);

   wfp=fopen(snnsFilenames.validationPatternsFilename,"w");

   write_patterns_into_file(wfp,
                            validationStartPointer,
                            numberOfValidationPatterns,
                            numberOfInputUnits,
                            numberOfOutputUnits);
   fclose(wfp);

   /* Step : Write a batchman script that does the following:

      1. Set parameters for various functions (learning, update etc.).
      2. Train the neural network.
      3. Validate the quality of the network by testing it against the validation patterns.
      4. Write a file in the format: actual angle fraction    predicted angle fraction
   */

   write_snns_commands(&snnsFilenames,&snnsParameters);

   /* Step : Execute the batchman script and calculate the Pearsons coefficient between the values
	      predicted by the neural network and the actual values.
   */

   execute_script_and_calculate_correlation_coefficient(validationStartPointer,&snnsFilenames);

   /* Step : Release memory allocated to the training and validation patterns. Also, remove
	      temporary files created.
   */

   free_list(); 
   remove_temporary_files();

   return 1;

} /* End of program. */
