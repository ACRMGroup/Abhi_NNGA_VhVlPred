# include "snns.h"
# define TEMP_PATH "/tmp/SNNS"

/* --------------------------- INTRODUCTION ---------------------------------


 This program performs the task of training, validating, and printing the correlation
 coefficient of the validation data to predicted data. The essential 

1. Default values:
		    -> init function - Randomize weights - Randomize_Weights (parameter to setInitFunc).
		    -> learning function - Resilient propogation - Rprop (parameter to setLearnFunc).
		    -> update function - Topological order - Topological_Order (parameter to setUpdateFunc).
		    -> pruning function - Magnitude based pruning - MagPruning (parameter to setPruningFunc).
		    -> shuffling mode - TRUE - ( setShuffle(TRUE) ).

2. Use default weights and parameters:

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


/* ------------------------ GLOBAL VARIABLES -------------------------------

static struct snns_filenames
{
   // Parameters pertaining to filenames 

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
   // Parameters for training the network

   int numberOfCycles;

   double sseThreshold;

   // Functions and associated parameters

   // setInitFunc(Randomize_Weights, 1.0, -1.0)

   char initFunction[20];
   double initFunctionParam1,initFunctionParam2;

   // setLearnFunc(Rprop,0.2,0,1e-6)

   char learningFunction[20];
   double learningFunctionParam1,
	 learningFunctionParam2;
   char learningFunctionParam3[8];

   // setUpdateFunc(Topological_Order)

   char updateFunction[20];

   // setPruningFunc("MagPruning", "Rprop", 15.0, 3.5, FALSE, 500, 90,1e6, 1.0)

   char pruningFunction[20];
   char pruningFunctionParam1[20];
   double pruningFunctionParam2,pruningFunctionParam3;
   char pruningFunctionParam4[8];
   double pruningFunctionParam5,pruningFunctionParam6;
   char pruningFunctionParam7[8];
   double pruningFunctionParam8;

   // setShuffle(TRUE)

   char shuffling[8];

   // Priority flag for over riding number of cycles or vice versa or using a combination
   // of user specified inputs for both

   int ssePriority;
   int cyclesPriority;

};

*/


/* Declare elements of type "snns_filenames" and "snns_parameters". */

static struct snns_filenames snnsFilenames;

static struct snns_parameters snnsParameters;

static struct pattern_variables patternVariables;

static int userNumberOfCycles = -1;

static double userSSEThreshold = -1;

static char snnsParametersConfigFilename[100];

static char *abhiDataDir=NULL;

static int numberOfInputUnits=0,
           numberOfOutputUnits=1;

static char intermediateFilename[100];

node *trainingStartPointer=NULL,
            *validationStartPointer=NULL;

static int currentProcessID=-1;

FILE *wfp=NULL;


/* ----------- END OF GLOBAL VARIABLE DECLARATION SECTION ------------------ */


/* --------------- SUB - ROUTINES DECLARATION SECTION ----------------------

   Functions defined in the library:
   ---------------------------------

   int read_intermediate_file(char *intermediateFilename,
                              struct pattern_variables *patternVariables,
                              node *trainingStartPointer,node *validationStartPointer);

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

   int write_snns_commands(struct snns_filenames *snnsFilenames,
                           struct snns_parameters *snnsParameters);

   BOOL execute_script_and_calculate_correlation_coefficient(node *validationFirst,
                                                             struct snns_filenames *snnsFilenames,
                                                             struct pattern_variables *patternVariables);

 
   BOOL initialise_snns_parameters(char *configFilename,struct snns_parameters *snnsParameters);
*/


void assign_snns_filenames(struct snns_filenames *snnsFilenames,int currentProcessID);

void free_list();

void remove_temporary_files();

void Usage();

BOOL parse_command_line_parameters(int numberOfParam,char **param);



/* ------------ END OF SUB - ROUTINES DECLARATION SECTION ------------------ */



/* ---------------- SUB - ROUTINES DEFINITION SECTION ---------------------- */



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

   unlink(snnsFilenames.snnsTrainingScriptFilename);
   unlink(snnsFilenames.snnsValidationScriptFilename);
   unlink(snnsFilenames.trainingPatternsFilename);
   unlink(snnsFilenames.validationPatternsFilename);
   unlink(snnsFilenames.trainedNetworkFilename);
   unlink(snnsFilenames.predictedValuesOutputFilename);
   unlink(snnsFilenames.inputOutputComparisonsFilename);

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
	 strcpy(patternVariables.trainingSetBoundaryString,param[i+1]);
	 i+=2;
	 parameterCount++;
	 continue;
      }
      else
      if(! strcasecmp(param[i],"-out") )
      {
	 strcpy(snnsFilenames.coefficientOutputFilename,param[i+1]);
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
   /* Step 1: Check if required number of inputs have been supplied by user. If not
	      print usage scheme of the program and exit.
   */

   if(argc < 6)
   {
      Usage(argv);
      return 0;
   }

   /* Step 2: Parse the command line parameters. If any problems of mismatch are
	      encountered, print usage scheme of program and exit.

	      Before entering the routine, make NULL, all parameters that are
	      optional in the program.
   */

   snnsParametersConfigFilename[0]=0;
   snnsFilenames.coefficientOutputFilename[0]=0;
   userSSEThreshold=-1;
   userNumberOfCycles=-1;

   if(! parse_command_line_parameters(argc,argv) )
   {
      Usage(argv);
      return 0;
   }

   numberOfInputUnits=read_intermediate_file(intermediateFilename,
					     &patternVariables,
					     &trainingStartPointer,&validationStartPointer);

   if(numberOfInputUnits < 0)
   {
      printf("Exiting program.\n\n");

      return 0;
   }

   currentProcessID=getpid();

   /* Step 3: Set parameters for the SNNS file names and SNNS training parameters */

   assign_snns_filenames(&snnsFilenames,currentProcessID);

   if(snnsParametersConfigFilename[0] == 0)
   {
      abhiDataDir=getenv("ABHIDATADIR");
      sprintf(snnsParametersConfigFilename,"%s/ABNUM/snns_training_parameters.txt",abhiDataDir);
   }

   if(! initialise_snns_parameters(snnsParametersConfigFilename,&snnsParameters) )
   {
      printf("Exiting program.\n\n");

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

   /* Step 4: Write the training and validation set patterns into respective files */

   wfp=fopen(snnsFilenames.trainingPatternsFilename,"w");

   write_patterns_into_file(wfp,
                            trainingStartPointer,
                            patternVariables.numberOfTrainingPatterns,
                            numberOfInputUnits,
                            numberOfOutputUnits);
   fclose(wfp);

   wfp=fopen(snnsFilenames.validationPatternsFilename,"w");

   write_patterns_into_file(wfp,
                            validationStartPointer,
                            patternVariables.numberOfValidationPatterns,
                            numberOfInputUnits,
                            numberOfOutputUnits);
   fclose(wfp);

   /* Step 5: Write a batchman script that does the following:

      1. Set parameters for various functions (learning, update etc.).
      2. Train the neural network.
      3. Validate the quality of the network by testing it against the validation patterns.
      4. Write a file in the format: actual angle fraction    predicted angle fraction
   */

   if(! write_snns_commands(&snnsFilenames,&snnsParameters) )
   {
      printf("Exiting program.\n\n");

      free_list();
      remove_temporary_files();

      return 0;
   }

   /* Step 6: Execute the batchman script and calculate the Pearsons coefficient between the values
	      predicted by the neural network and the actual values.
   */

   if(! execute_script_and_calculate_correlation_coefficient(validationStartPointer,&snnsFilenames,&patternVariables) )
   {
      printf("Exiting program.\n\n");

      free_list();
      remove_temporary_files();
 
      return 0;
   }

   /* Step 7: Release memory allocated to the training and validation patterns. Also, remove
	      temporary files created.
   */

   free_list(); 
   /* remove_temporary_files(); */

   return 1;

} /* End of program. */
