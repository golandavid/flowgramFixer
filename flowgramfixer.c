// flowgramfixer.c-- convert flowgram incorporation sequences into nucleotide
//                   sequences
//
// References:
//
//	[1]	Golan and Medvedev, "Using State Machines to Model the Ion Torrent 
//		Sequencing Process and Improve Read Error-Rates" (work in progress?)

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
//#define NUM_THREADS   10

#define negInf (-DBL_MAX)

#define CONST 1 
#define INTERCEPT 2 
#define TREND 3 
#define TREND_ONLY 4 
#define GREEDY 5 


#define NORMAL 0
#define EXP   1

/*
#ifndef M_PI
#define M_PI 3.14159
//#define M_PI 3.141592653589793
//#define M_PI 3.141592654
#endif
*/
//----------
//
// debugging switches--
//
//----------

//#define debugObssLength
//#define debugSeqLength

//----------
//
// global data and types--
//
//----------

// states for the state machine model
//
// A state is (conceptually) a set of the four nucleotides;  see the section
// "operations on single states" for more detail.
//
// We define types for vectors indexed by the 16 possible states.  These are
// used for rows in the viterbi algorithm's DP matrix.

#define NUM_STATES 16

typedef char sstate[5];		// a string representation of the form "AC-T" or ""
typedef int    introw[NUM_STATES];
typedef double realrow[NUM_STATES];
typedef struct struct_it{
 double intercept;
 double trend;

} struct_it;

// Struct inludes preprocessed likelihood matrix with and without incorporation
typedef struct prep_matrix{
	int*** whichMax_matrix;
	double*** lik_with_incorp_matrix;
	double** lik_without_incorp_matrix;
	}prep_matrix;

// Struct includes all the values passed to the threads.
typedef struct thread_data{
	char* chunk_in;
	char* chunk_out;
	char* washSequence;
	int mode;
	int dist;
	double sd_trend;
	double sd_start;
	double step_size;
	int maxlen;
	char* ml_flag;
	prep_matrix prep;
	}thread_data;


// state for reading text files line-by-line; see read_line()
typedef struct filestate {
	char*	filename;
	int		lineNum;
	int		missingEol;
	} filestate;

// other globals

double TRANS_PROBS[NUM_STATES][NUM_STATES][4];
char* defaultWashSequence = "tacgtacgtctgagcatcgatcgatgtacagc";

int dist ; 
int mode ; 
double sd_start;
double sd_trend ; 

int UNDERFLOW_CONST=2;

double gcContent;

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

static void	usage			(char* programName);
prep_matrix 	generate_lh_table	( double intercept, double trend, int maxlen, int dist);
void		get_ml_estimates	(char * filename, double step_size, char* washSequence,  double *sample_intercept,double *sample_trend);
void*		thread_call		( void *thread_arg);
char* 		concat			(char *s1, char *s2);
void 		split_sff 		(int parts, char* filename, char* chunk_names[], double percent, int* maxLen);
unsigned long 	count_lines		(char *filename);
void 		join_seq		(char* out_file, int n);
void 		join_lik		(char* out_file, int n);
void 		delete_chunks		(char* out_file, int n);
void 		delete_intermediates	(char* out_file, int n);
void 		fill_const		(double* params, int paramsLen, double fill);
void 		fill_trend		(double* params, int paramsLen, double intercept, double trend);
struct_it	return_optimal_trend_intercept_greedy (double* obss, int obssLen, double* params, char* washSequence, double int_start, double trend_start, double int_step, double trend_step);
void 		find_optimal_const_sigma (double* obss,int obssLen,double* params, char* washSequence, double start, double end, int steps);
double 		find_optimal_trend_sigma (double* obss,int obssLen,double* params, char* washSequence, double int_start, double int_end, double trend_start, double trend_end, int int_steps, int trend_steps);
double 		find_optimal_trend_only	 (double* obss,int obssLen,double* params, char* washSequence, double intercept, double trend_start, double trend_end, int trend_steps);
double 		find_optimal_trend_greedy(double* obss,int obssLen,double* params, char* washSequence, double int_start, double trend_start, double int_step, double trend_step);
static int*	viterbi			(double* obss, int obssLen, double* params, char* washSequence, int* seqLen, int dist);
static int*  	viterbi_prepros		(double* obss, int obssLen, double* params, char* washSequence, int* seqLen, prep_matrix prep, char* ml_tag, int dist);
static double	forward			(double* obss, int obssLen, double* params, char* washSequence);
static void	create_transition_probs	(double gcContent);
static int	index_of_max		(realrow vals);
static void	one_step		(double obs, char wash, int washIndex, const realrow oldScores, realrow newScores, introw newBest, introw newBestInc, double param, int dist);
static void	one_step_prepros	(double obs, char wash, int washIndex, const realrow oldScores, realrow newScores, introw newBest, introw newBestInc, double param, prep_matrix prep, char* ml_tag, int dist);
static void	one_step_forward	(double obs, char wash, int washIndex, const realrow oldProbs, realrow newProbs, double params);
static double	obs_given_no_incorp	(double obs, double param,int dist);
static double	obs_given_incorp	(double obs, char wash, int* whichMax, double param,int dist);
static double	obs_given_incorp_forward(double obs, char wash, int washIndex, double param, int dist);
static double	dnorm			(double x, double mean, double sd);
static double	ddexp			(double x, double mean, double sd);
static double	dgeom			(int x, double prob);
static void	sstate_from_num		(sstate str, int statenum);
static int	sstate_to_num		(sstate str);
static int	sstate_get		(sstate str, int wash);
static void	sstate_see_char		(sstate str, char c);
static void	sstate_see_num		(sstate str, int x);
static void	sstate_skip_char	(sstate str, char c);
static void	sstate_skip_num		(sstate str, int x);
static int	nt2num			(char c);
static FILE*	open_file		(char* basename, const char* extension, const char* mode);
static void	init_read_line		(filestate* fState, char* filename);
static int	read_line		(FILE* f, filestate* fState, char* line, size_t lineSize);
static double*	get_observations	(FILE* f, filestate* fState, int* obssLen);
static void	zero_trim_observations	(double* obss, int* obssLen);
static void	write_flow_as_nts	(FILE* f, int* seq, int seqLen, char* washSequence);
static double	string_to_double	(const char* s);
static char*	skip_whitespace		(char* s);
static char*	skip_darkspace		(char* s);
void 		print_time		();

//----------
//
// main program--
//
//----------


// function describes the usage
static void usage
   (char* programName)
	{
	fprintf (stderr, "usage: %s -i incorporation_file -o output_file -d dist -m mode \n", programName);
	fprintf (stderr, "\n");
	fprintf (stderr, " -i \t incorporation_file is space delimited file of incorporation sequences\n");
	fprintf (stderr, " -o \t output_file is the prefix for all output files\n");
	fprintf (stderr, " -d \t exp/normal distribution to be used in the likelihood, (default value is normal)\n");
	fprintf (stderr, " -m \t greedy/trend/const , algorithm being used(greedy is greedy search, trend is grid search, const is user defined), default is greedy\n");
	fprintf (stderr, " -s \t additional parameter for greedy algorithm (step size - default is 1, larger step size runs faster) \n");
	fprintf (stderr, " -p \t additional parameter for greedy algorithm, percentage of sequences used to generate the intercept and trend (default value is 1)\n");
	fprintf (stderr, " -x \t set when mode is const, the  intercept for the noise model of all the reads \n");
	fprintf (stderr, " -y \t set when mode is const, the  trend for the noise model of all the reads (default is 0)\n");	
	fprintf (stderr, " -t \t number of threads used to run the program (default is 1)\n");
	fprintf (stderr, " -l \t generate preprocessed maximum likelihood values (default is YES)\n");
	fprintf (stderr, " -w \t define the washcycle in lower case, default value is 'tacgtacgtctgagcatcgatcgatgtacagc' \n");
	fprintf (stderr, " -v \t YES/NO/ADVANCE; YES prints status info, ADVANCE prints status info and keeps all intermediate files, NO does not print any info, useful to debug (default is NO)\n");
	fprintf (stderr, "\n");
	exit (EXIT_FAILURE);
	}


//--- for debugObssLength ---

//#define debugObssLength

#ifndef debugObssLength
#define debugObssLength_1 ;
#define debugObssLength_2 ;
#endif // not debugObssLength

#ifdef debugObssLength

#define debugObssLength_1                                                     \
	{                                                                         \
	int ix;                                                                   \
	fprintf (stderr, "read %d obss (%d):", readNumber, obssLen);              \
	for (ix=0 ; ix<obssLen ; ix++)                                            \
		fprintf (stderr, " %0.2f",obss[ix]);                                  \
	fprintf (stderr, "\n");                                                   \
	}

#define debugObssLength_2                                                     \
	{                                                                         \
	int ix;                                                                   \
	fprintf (stderr, "trimmed obss (%d):", obssLen);                          \
	for (ix=0 ; ix<obssLen ; ix++)                                            \
		fprintf (stderr, " %0.2f",obss[ix]);                                  \
	fprintf (stderr, "\n");                                                   \
	}

#endif // debugObssLength


//--- for debugSeqLength ---

//#define debugSeqLength

#ifndef debugSeqLength
#define debugSeqLength_1 ;
#endif // not debugSeqLength

#ifdef debugSeqLength

#define debugSeqLength_1                                                      \
	{                                                                         \
	int ix;                                                                   \
	fprintf (stderr, "seq (%d):", seqLen);                                    \
	for (ix=0 ; ix<seqLen ; ix++)                                             \
		fprintf (stderr, " %0d",seq[ix]);                                     \
	fprintf (stderr, "\n");                                                   \
	}

#endif // debugSeqLength

//--- main ---
int  main (int argc, char **argv)
{
  
  
  //initialize and read the command line parameters
  char *inputflag = NULL;
  char *outputflag = NULL;
  char *dist_flag = "normal";
  char *mode_flag = "greedy";
  char *percent_flag = "1";
  char *threadflag= "1";
  char *stepsize_flag="1";
  char *intercept_flag=NULL;
  char *trend_flag=NULL;
  char *ml_matrix_flag="YES";
  char *washcycle_flag=defaultWashSequence;
  char *verbose_flag="NO";
  int index;
  int c;
  opterr = 0;
  
  // reading all the arguments frrom the command line
  while ((c = getopt (argc, argv, "i:o:d:m:p:t:x:y:s:l:w:v:")) != -1){
    switch (c)
    {
      case 'i':
        inputflag = optarg;
        break;
      case 'o':
        outputflag = optarg;
        break;
      case 'd':
        dist_flag = optarg;
        break;
      case 'm':
        mode_flag = optarg;
        break;
      case 'p':
        percent_flag = optarg;
        break;
      case 't':
        threadflag = optarg;
        break;
      case 'x':
        intercept_flag = optarg;
        break;
      case 'y':
        trend_flag = optarg;
        break;
      case 's':
        stepsize_flag = optarg;
        break;
      case 'l':
        ml_matrix_flag = optarg;
        break;
      case 'w':
        washcycle_flag = optarg;
        break;
      case 'v':
        verbose_flag = optarg;
        break;
      case '?':
        if (optopt == 'i') {
            fprintf (stderr, "Option -%c missing an argument, please enter input file name.\n", optopt);
            abort ();
        } else if (optopt == 'o') {
            fprintf (stderr, "Option -%c missing an argument, please enter the outfile name.\n", optopt);
            abort ();
        } else if (optopt == 'd') {
            fprintf (stderr, "Option -%c missing an argument, please set value as normal/exp.\n", optopt);
            abort ();
        } else if (optopt == 'm') {
            fprintf (stderr, "Option -%c missing an argument, please set value as greedy/trend/const.\n", optopt);
            abort ();
        } else if (optopt == 'p') {
            fprintf (stderr, "Option -%c missing an argument, can ommit the flag to set the default value.\n", optopt);
            abort ();
        } else if (optopt == 't') {
            fprintf (stderr, "Option -%c missing an argument, can ommit the flag to set the default value.\n", optopt);
            abort ();
        } else if (optopt == 's') {
            fprintf (stderr, "Option -%c missing an argument, can ommit the flag to set the default value.\n", optopt);
            abort ();
        } else if (optopt == 'w') {
            fprintf (stderr, "Option -%c missing an argument, can ommit the flag to set the default value.\n", optopt);
            abort ();
        } else if (optopt == 'x') {
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            abort ();
        } else if (optopt == 'y') {
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            abort ();
	} else if (optopt == 'l') {
            fprintf (stderr, "Option -%c missing an argument, requires YES/NO.\n", optopt);
            abort ();
        } else if (optopt == 'v') {
            fprintf (stderr, "Option -%c missing an argument, requires YES/NO.\n", optopt);
            abort ();
        } else if (isprint (optopt)) {
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            abort ();
        } else {
            fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
            default:
            abort ();
        }
    }
  }
  for (index = optind; index < argc; index++){
    printf ("Non-option argument %s\n", argv[index]);
  }       
  
  //convert the command line parameter variables to required datatypes
  char*	 washSequence;
  char*	 allObssFilename;
  char*	 outputTemplate;
  int	 readNumber;
  char*  ml_tag;
  int  verb_tag;
  double step_size;  
	
  if (inputflag == NULL) {
    printf ("ERROR: -i flag is missing\n");
    usage (argv[0]);
    abort ();
  }
  if (outputflag == NULL) {
    printf ("ERROR: -o flag is missing\n");
    usage (argv[0]);
    abort ();
  }
	
  
  if( access( inputflag, F_OK ) != -1 ) {
    allObssFilename = inputflag;
    // file exists
  } else {
    // file doesn't exist
    printf("ERROR: Input file not found\n");
    usage (argv[0]);
    abort ();
  }
  outputTemplate  = outputflag;
  ////parse dist argument 
  if (strcmp(dist_flag,"normal")==0) {
      dist=NORMAL; 
  } else if (strcmp(dist_flag,"exp")==0) {
      dist=EXP; 
  } else {
      printf ("ERROR: -d flag invalid value\n");
      usage (argv[0]);
      abort ();
  }
  
  //// parse param arg, which gives us the mode
  if(strcmp(mode_flag,"intercept")==0){
      mode=INTERCEPT; 
  } else if (strcmp (mode_flag,"trend") == 0){
      mode=TREND; 
  } else if (strcmp (mode_flag,"trendOnly") == 0) {
      mode=TREND_ONLY;
      sd_start=string_to_double (intercept_flag);
  } else if (strcmp (mode_flag,"greedy") == 0 ) {
      mode=GREEDY ; 
      step_size = string_to_double (stepsize_flag);
  } else if (strcmp (mode_flag,"const") == 0) {
      mode=CONST; 
      if(trend_flag!=NULL){
        sd_start  = string_to_double (intercept_flag);
        if (intercept_flag!=NULL) {
	  sd_trend  = string_to_double (trend_flag);
        } else {
	  sd_trend = 0 ; 
        }
      } else {
        printf("ERROR: Missing parameter's -x\n");
        usage(argv[0]);
        abort ();
      }  
  } else {
      printf("ERROR: -m flag invalid value\n");
      usage(argv[0]);
      abort ();
  }
  //// parse the maximum likelihood matrix flag   
  if (strcasecmp (ml_matrix_flag,"YES") == 0) {
      ml_tag="YES"; 
  } else if (strcasecmp (ml_matrix_flag,"NO") == 0) {
      ml_tag="NO"; 
  }else {
      printf("ERROR: -l flag invalid value\n");
      usage(argv[0]);
      abort ();
  }
 //// parse the verbose_flag
    if (strcasecmp (verbose_flag,"YES") == 0) {
      verb_tag=1; 
  } else if (strcasecmp (verbose_flag,"NO") == 0) {
      verb_tag=0; 
  } else if (strcasecmp (verbose_flag,"ADVANCE") == 0) {
      verb_tag=2; 
  }else {
      printf("ERROR: -v flag invalid value\n");
      usage(argv[0]);
      abort ();
  }

  //// print the mode and distribution used.	
  if (verb_tag>=1) {
   print_time();
   printf("Using the %s distribution (%d) with mode %s (%d)\n",dist_flag,dist,mode_flag,mode);
  }
  //// can add gcContenet and washSequence as parameters in fututre
  gcContent = 0.5;
  washSequence = washcycle_flag;
  //// parse the number of threads and initialize the threads
  int NUM_THREADS=atoi(threadflag);
  thread_data td[NUM_THREADS];

  
  //// These variables are used in generation of preprocessed ML matrix.
  prep_matrix prep;
  double sample_intercept;
  double sample_trend;
  //// parse the percentage value to generate subsample
  double percent_seq=atof(percent_flag);
  //Splitting the input flowgram into smaller files and create a subsample of the input dataset
  // the split files created have the same input file name with attached "_0NN" where NN is the number of threads 01 to 10.
  // the subsample file has the same input file name followed by subsample
  char* sffchunks[NUM_THREADS];
  int maxLen=0;
  if (verb_tag>=1) {
    printf("Splitting Input file...\n ");
  }
  
  split_sff(NUM_THREADS,allObssFilename,sffchunks, percent_seq,&maxLen);
  if (verb_tag>=1) {
   printf("max length of flow= %d\n",maxLen);
   print_time();
  }
  //Obtain the mean intercept and trend from subsample of the input flowgram.
  // implemented only for GREEDY mode.

  if (mode==GREEDY) {
    if (verb_tag>=1) {
     printf("Generating Likelihood Matrix....\n ");
     printf ("Averaged from %f%% of the data: \n",percent_seq);
    } 
    get_ml_estimates(allObssFilename,step_size,washSequence,&sample_intercept,&sample_trend);
    if (verb_tag>=1) {  print_time(); }
    sd_start=sample_intercept;
    sd_trend=sample_trend;
   //// the program uses the trend and intercept from subsample and process the whole input file.
   //// mode is changed to CONST since trend and intercept are known. 
    mode=CONST;
  }
  
  if (mode==CONST) {
    if (verb_tag>=1) { printf("intercept\t%.11f\ntrend\t%.11f\n",sd_start,sd_trend); }
   //// The maximum likelihood values for all possible observations 0.01 to 12.99 were caliculated at different positions in the flowgram.
    prep=generate_lh_table(sd_start,sd_trend,maxLen,dist);  
  }

  int rc, id, idx;
  char chunk_f[NUM_THREADS][100], chunk_o[NUM_THREADS][100];
  
  //initialize the threads and make them joinable
  pthread_t threads[NUM_THREADS];
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (verb_tag>=1) { printf("Processing using threads.....\n "); }	
  //Perform base calling using multiple threads
  for (id=0; id<NUM_THREADS; id++) {
    char * washSeq= washSequence;	     
    //generating the input and output files names passed to thread.
    sprintf(chunk_f[id], "%s_%03d", allObssFilename,id);
    sprintf(chunk_o[id],"%s_%03d", outputTemplate,id);
    //assign values to struct  
    td[id].chunk_in = chunk_f[id];
    td[id].washSequence = washSeq;
    td[id].sd_trend = sd_trend;
    td[id].mode = mode;
    td[id].dist = dist;		
    td[id].sd_start = sd_start;
    td[id].chunk_out= chunk_o[id];
    td[id].step_size=step_size;
    td[id].maxlen=maxLen;
    td[id].ml_flag=ml_tag;
    td[id].prep=prep;

    //Creating the threads and calling the function
    rc = pthread_create(&threads[id], &attr, thread_call, (void *) &td[id] );
    if (rc) {
      printf("Error:unable to create thread\n");
      exit(-1);
    }
  }
  
  //Joining the threads
  pthread_attr_destroy(&attr);
  for (idx=0; idx<NUM_THREADS; idx++) {
    rc = pthread_join(threads[idx], NULL);
    if (rc) {
      printf("Error:unable to join\n");
      exit(-1);
    }
  }      

  //Joining the output files generated by each individual thread		
  
  join_seq(outputTemplate,NUM_THREADS);   
  join_lik(outputTemplate,NUM_THREADS);
  if (verb_tag>=1) {  print_time(); }
  //Deleting the intermediate files
  if (verb_tag<2){
    delete_chunks(allObssFilename,NUM_THREADS);
    delete_intermediates(outputTemplate,NUM_THREADS);
  }
  pthread_exit(NULL);
 
  return EXIT_SUCCESS;
}


//From the one percent subsample generated using split_sff(), the intercept and trend are generated using greedy approach(return_optimal_trend_intercept_greedy()) for each flow
//and the mean of all the intercepts and trend is updated. 
//
void get_ml_estimates (char * filename, double step_size, char* washSequence, double *sample_intercept,double *sample_trend) {

  char fname[100];
  double sum_i=0, sum_t=0;
  ////reads the subsample generated in split_sff()
  sprintf(fname, "%s_subsample", filename);
  FILE *subset = fopen(fname, "rt");
  create_transition_probs (gcContent);
  int	readNumber = 0; 
  char	line[1000*6+1];
  //// reads each line and caliculates the intercept and trend for each flow
  while (fgets (line, sizeof(line), subset) != NULL) {
    double* obss;
    double*    params;
    int		obssLen;
    int		seqLen;
    double lik=0 ;
    //double intercept_return;
    //double trend_return; 
    char *saveptr1;
    const char s[2] = " ";
    char* temp_num;
    int e_num = 0;
    struct_it pass;
    obss = malloc(sizeof(double)*sizeof(line));
    temp_num = strtok_r(line,s,&saveptr1);//split string
    ////splits the flowgram by white sapce and reads each value one after the other
    while (temp_num != NULL) {
      double tmp_d=atof(temp_num);
      obss[e_num] = tmp_d;
      temp_num = strtok_r(NULL,s,&saveptr1);
      e_num++;
    }
    obssLen=e_num;
    if (obss == NULL) break;
    
    readNumber++;
    debugObssLength_1;
    zero_trim_observations (obss, &obssLen);
    debugObssLength_2;
    params = malloc(obssLen * sizeof(double));
    pass=return_optimal_trend_intercept_greedy(obss, obssLen, params, washSequence,0.03,0.0002,0.001*step_size,0.00001*step_size);
//find_optimal_trend_greedy(obss, obssLen, params, washSequence,0.03,0.0002,0.001*step_size,0.00001*step_size);
    sum_i+=pass.intercept;
    sum_t+=pass.trend;
    free (params);
    free (obss);
 }
   
 *sample_intercept=(sum_i/readNumber);
 *sample_trend=(sum_t/readNumber);
 
}

//populates the matrix with maximum likelihood values for observations values 0 to 13, gcContent=0.5( all 4 nuclotides have same probability of occurance)
//and for all the positions in the flowgram..   
prep_matrix generate_lh_table(double intercept, double trend, int maxlen, int dist)
{
  prep_matrix prep;
  int i,j,k,in,jn;  
  //the order of nucleotides matches the order in nt2num(), in future make sure the order ACGT is same as in nt2num function.
  char wash[4]={'A','C','G','T'};
  prep.lik_with_incorp_matrix=(double***) malloc(1300*sizeof(double**));
  prep.whichMax_matrix=(int***) malloc(1300*sizeof(int**));
  prep.lik_without_incorp_matrix=(double**) malloc(1300*sizeof(double*));
  for (in=0;in<1300;in++) {
    prep.lik_with_incorp_matrix[in]=(double**)malloc(4*sizeof(double*));
    prep.whichMax_matrix[in]=(int**)malloc(4*sizeof(int*));
    prep.lik_without_incorp_matrix[in]=(double*)malloc(maxlen*sizeof(double));
    for(jn=0;jn<4;jn++){
      prep.lik_with_incorp_matrix[in][jn]=(double*)malloc(maxlen*sizeof(double));
      prep.whichMax_matrix[in][jn]=(int*)malloc(maxlen*sizeof(int));
    }
  }
  for(i=0;i<1300;i++){
    for(j=0;j<4;j++){
      for(k=0;k<maxlen;k++){
        double obs=(double)i/100;
        double sd = intercept + ((double)k)*trend ; 
        int temp_max=0;
        double temp=obs_given_incorp(obs,wash[j],&temp_max,sd,dist);
        prep.lik_with_incorp_matrix[i][j][k]=temp;
        prep.whichMax_matrix[i][j][k]=temp_max;
      }
    }
  }
 
  int l,m;
  for(l=0;l<1300;l++){
    for(m=0;m<maxlen;m++){
      double obs_val=(double)l/100;
      double sd_val = intercept + ((double)m)*trend ; 
      double temp_val=obs_given_no_incorp(obs_val,sd_val,dist);
      prep.lik_without_incorp_matrix[l][m]=temp_val;
    }
  }
	
  return prep;
} 

// This is the function which is called in each thread
// It process each chunk independently writes the resulting output to output_0NN.lik and outpt_0NN.fa files
// 
void *thread_call( void *targ)
{
  //the GC content is set to 0.5 by default
  double gcContent = 0.5;
  create_transition_probs (gcContent);
  //the struct thread_data which includes the parameters sent to each thread is initiated. 
  thread_data* arg= (thread_data*) targ;
  filestate	allObssState;
  int			readNumber;
  int			i,j;
  char	line[1000*6+1];
  
  //reads the input file and initiates the file handel for output files.
  FILE *allObssF,*likF,*seqF;
  allObssF = open_file (arg->chunk_in, NULL, "rt");
  init_read_line (&allObssState, arg->chunk_in);
  likF    = open_file (arg->chunk_out, ".lik", "wt");
  seqF     = open_file (arg->chunk_out, ".fa",  "wt");
  readNumber = 0; 
  //for each line in the input file perform base call based on the mode specified.
  while (fgets (line, sizeof(line), allObssF) != NULL) {
    double* obss;
    double* params;
    int*    seq;
    int	    obssLen;
    int	    seqLen;
    double  lik ;
    char    *saveptr1;
    const char s[2] = " ";
    char*   temp_num;
    lik=0;
    obss = malloc(sizeof(double)*sizeof(line));
    temp_num = strtok_r(line,s,&saveptr1);//split string
    int     e_num = 0;
    while(temp_num != NULL) {
      double tmp_d=atof(temp_num);
      obss[e_num] = tmp_d;
      temp_num = strtok_r(NULL,s,&saveptr1);
      e_num++;
    }

    obssLen=e_num;
    if (obss == NULL) break;
    readNumber++;
    debugObssLength_1;
    zero_trim_observations (obss, &obssLen);
    debugObssLength_2;
    params = malloc(obssLen * sizeof(double));
    //// checks for the input mode and then calls the viterbi function accordingly
    if (arg->mode==CONST) {
      fill_trend(params,obssLen,arg->sd_start,arg->sd_trend) ; 
      ////viterbi_process is identical to viterbi in functionality. The only difference is it uses the preprocessesd ml values instead of processing them at each step 
      seq = viterbi_prepros (obss, obssLen, params, arg->washSequence, &seqLen, arg->prep, arg->ml_flag, arg->dist);
    } else if (arg->mode==INTERCEPT) {
       find_optimal_const_sigma(obss, obssLen, params, arg->washSequence,0,0.2,10);
       seq = viterbi (obss, obssLen, params, arg->washSequence, &seqLen, arg->dist);
    } else if (arg->mode==TREND_ONLY){
       lik=find_optimal_trend_only(obss, obssLen, params, arg->washSequence,arg->sd_start,0,0.0005,10);
       seq = viterbi (obss, obssLen, params, arg->washSequence, &seqLen, arg->dist);
    } else if (arg->mode==GREEDY) {
       lik=find_optimal_trend_greedy(obss, obssLen, params, arg->washSequence,0.03,0.0002,0.001*arg->step_size,0.00001*arg->step_size);
       seq = viterbi (obss, obssLen, params, arg->washSequence, &seqLen, arg->dist);
    } else { // mode is trend 
       lik=find_optimal_trend_sigma(obss, obssLen, params, arg->washSequence,0,0.1,0,0.0005,10,20);
       seq = viterbi (obss, obssLen, params, arg->washSequence, &seqLen, arg->dist);
    }
    
    fprintf(likF,"%f\t%f\t%f\n",params[0],params[1]-params[0],lik);
    debugSeqLength_1
    fprintf(seqF,">File: %s Line: %d\n", arg->chunk_in,readNumber);
    write_flow_as_nts (seqF, seq, seqLen, arg->washSequence);
    free (params);
    free (obss);
    free (seq);
    fflush(likF);
    fflush(seqF);	
  }
  fclose (allObssF);
  fclose (likF);
  fclose (seqF);
  pthread_exit(0);
  return 0;
}

//Splits the input flowgram file into smaller chunks, the number of chunks is equal to the number of threads.
//Generates a subset file which includes 1% of the sequences taken from equal intervals in the input flogram.
void split_sff (int parts, char* filename, char* chunk_names[], double percent, int* maxLen )
{
  unsigned long line_count=count_lines(filename);
  double chunk_size= ceil (line_count/parts);
  double diff_count= line_count%parts;
  double lines_per_percent=ceil((line_count*percent)/100);
  double interval_line=ceil(line_count/lines_per_percent);
  int j,k,l;
  FILE *ss_infile=fopen(filename, "rt");
  char filename_tmp[100];
  sprintf(filename_tmp, "%s_subsample", filename);
  FILE *subfile=fopen(filename_tmp,"wt");
  int ss_size = 1024;
  int ss_pos;
  int ss_c;
  char *ss_buffer = (char *)malloc(ss_size);
  int check=0; 
  for(l=0;l<line_count;l++){
    ss_pos = 0;
    do { // read one line
      ss_c = fgetc(ss_infile);
      if (ss_c != EOF) {
        ss_buffer[ss_pos++] = (char)ss_c;
      }   
      if(ss_pos >= ss_size - 1) { // increase buffer length - leave room for 0
          ss_size *=2;	
          ss_buffer = (char*)realloc(ss_buffer, ss_size);
       }
    } while(ss_c != EOF && ss_c != '\n');

    check++;
    ss_buffer[ss_pos] = 0;
    if ((check+1)==interval_line) {
      fputs(ss_buffer,subfile);		   	
      check=0;
    }	 
  }
  fclose(ss_infile);
  fclose(subfile);
  
  FILE *insff=fopen(filename, "rt");
  if (insff == NULL) {
    exit(EXIT_FAILURE);	
  } 
  double no_lines;
  for ( j=0; j<parts; j++) {
    char filename_chunk[100];
    sprintf(filename_chunk, "%s_%03d", filename,j);
    size_t len = 0;
    chunk_names[j]=filename_chunk;
    FILE *sffchunk=fopen(filename_chunk, "wt");
    if (sffchunk == NULL) {
      exit(EXIT_FAILURE);	
    }   
    if (j==(parts-1)) {
      no_lines=chunk_size+diff_count;
    } else {
       no_lines=chunk_size;
    } 
    char * line = NULL;
    size_t lent = 0;
    ssize_t read;
    int pos, size = 1024;
    int c;
    char *buffer = (char *)malloc(size);
    for(k=0; k<no_lines; k++) {
      pos = 0;
      do { // read one line
        c = fgetc(insff);
	if (c != EOF) { buffer[pos++] = (char)c; }  
	if (pos >= size - 1) { // increase buffer length - leave room for 0
	  size *=2;	
	  buffer = (char*)realloc(buffer, size);
        }
      } while(c != EOF && c != '\n');
      buffer[pos] = 0;
      fputs(buffer,sffchunk);		   	
      ////keeps track of the lenght of each flow and gets the count for the longest flowgram
      int obssLen=max_obs_count(buffer);
      if (obssLen>=*maxLen) {
        *maxLen=obssLen;
      }
    }
    fclose(sffchunk);
   }
   fclose(insff);
}	

////gets the count of observations in each line
int max_obs_count (char * line) {
  char *saveptr1;
  const char s[2] = " ";
  char* temp_num;
  int e_num = 0;
  temp_num = strtok_r(line,s,&saveptr1);//split string
  ////splits the flowgram by white sapce and reads each value one after the other
  while (temp_num != NULL) {
    temp_num = strtok_r(NULL,s,&saveptr1);
    e_num++;
  }
  return e_num;
} 

// count the number of lines
unsigned long count_lines(char *filename) {
  FILE *fh = fopen(filename, "r");
  unsigned long line_count = 0;
  int c;
  c=getc(fh);
  while(c!= EOF) {
    if(c=='\n'){
      line_count++;
    }
    c=getc(fh);     
  }
  return line_count;
}


// combine all the .fa files from different threads into one large file.
void join_seq(char* out_file, int n)
{
  FILE *OUT;
  char out_name[100];
  sprintf(out_name,"%s.fa",out_file);
  OUT=fopen(out_name,"wt");
  if ( OUT == NULL ) {
    perror("Error ");
    exit(EXIT_FAILURE);
  } 
  int i;
  for(i=0; i<n; i++){
    char ch,filename_chunk[100];
    sprintf(filename_chunk, "%s_%03d.fa", out_file,i); // change ext for different files
    FILE * ft;
    ft = fopen(filename_chunk,"rt");
    if ( ft == NULL ) {
      perror("Error ");
      exit(EXIT_FAILURE);
    } 	
    while( ( ch = fgetc(ft) ) != EOF ) {
     fputc(ch,OUT);
    }
    fclose(ft);
  }
  fclose(OUT);
}

// combine all the .lik files from different threads into one large file.
void join_lik(char* out_file, int n)
{
  FILE *OUT;
  char out_name[100];
  sprintf(out_name,"%s.lik",out_file);
  OUT=fopen(out_name,"wt");
  if( OUT == NULL )
   {
      perror("Error ");
      exit(EXIT_FAILURE);
   } 
  int i;
  for(i=0; i<n; i++){
   char ch,filename_chunk[100];
   sprintf(filename_chunk, "%s_%03d.lik", out_file,i); // change ext for different files
   FILE * ft;
   ft = fopen(filename_chunk,"rt");
    
   if( ft == NULL )
   {
      perror("Error ");
      exit(EXIT_FAILURE);
   } 	
   while( ( ch = fgetc(ft) ) != EOF )
   {
     fputc(ch,OUT);
   }
   fclose(ft);
  }
  fclose(OUT);

}

//deletes the flowgram chunks created by split_sff() function.
void delete_chunks(char* out_file, int n)
{
  int i;
  for(i=0; i<n; i++){
    char seq_fn[100];
    sprintf(seq_fn,"%s_%03d",out_file,i);
    if( remove(seq_fn) != 0){
     printf("Unable to delete chunk files: %s\n", seq_fn);
    }
    
  }
  char seq_op[100];
  sprintf(seq_op,"%s_subsample",out_file);
    if( remove(seq_op) != 0){
     printf("Unable to delete subsample file: %s\n", seq_op);
    }
    

}

//deletes all the intermediate files generated by different threads.
void delete_intermediates(char* out_file, int n)
{
  int i;
  
  for(i=0; i<n; i++){
    char seq_fn[100];
    sprintf(seq_fn,"%s_%03d.fa",out_file,i);
    char lik_fn[100];
    sprintf(lik_fn,"%s_%03d.lik",out_file,i);
    if( remove(seq_fn) != 0){
     printf("Unable to delete .fa file: %s\n",seq_fn);
    }
    if(remove(lik_fn) != 0){
     printf("Unable to delete .lik file: %s\n",lik_fn);
    } 
  }

}

void fill_const
	(double* params,
	 int paramsLen,
	 double fill){
		int i ; 
		for(i = 0 ; i < paramsLen ; i++){
			params[i] = fill ; 
		}
	}

	
void fill_trend
	(double* params,
	 int paramsLen,
	 double start,
	 double trend)
	 {
		int i ; 
		
		for(i = 0 ; i < paramsLen ; i++){
			params[i] = start + ((double)i)*trend ; 
		}
	}

	
void find_optimal_const_sigma
	(double* obss,
	 int obssLen,
	 double* params, // this is the output 
	 char* washSequence,
	 double start,
	 double end,
	 int steps){
	
	double max_lik,maximizer,curr_sigma,lik ; 
	maximizer=start;
	max_lik=-1; 
	for(curr_sigma = start ; curr_sigma <= end ; curr_sigma += (end-start)/steps){
		// fill the params vector 
		fill_const(params, obssLen, curr_sigma); 

		lik = forward(obss,obssLen,params,washSequence);



		if(lik > max_lik){
			maximizer=curr_sigma; 
			max_lik = lik ; 
//			printf("%f\t%f",curr_sigma,lik);
//			printf("***");  
		}
//		printf("\n");
	}
	// fill with the maximizer 
	fill_const(params, obssLen, maximizer); 
	return;
}	

	
double find_optimal_trend_sigma
	(double* obss,
	 int obssLen,
	 double* params, // this is the output 
	 char* washSequence,
	 double int_start,
	 double int_end,
	 double trend_start,
	 double trend_end,
	 int int_steps,
	 int trend_steps){
	
	double max_lik,curr_sigma,lik ;
	double step1,step2;
	double intercept,trend; 
	double maximizers[2] ; 
	step1 = (int_end-int_start)/int_steps; 
	step2 = (trend_end-trend_start)/trend_steps ; 
	
	maximizers[0]=-1;
	maximizers[1]=-1;
	max_lik=-1; 
	for(intercept = int_start ; intercept <= int_end ; intercept += step1){
		for(trend = trend_start ; trend <= trend_end ; trend += step2){
		// fill the params vector 
		fill_trend(params, obssLen, intercept,trend); 
		lik = forward(obss,obssLen,params,washSequence);

			printf("%f-%f\t%f",intercept,trend,lik);
			if(lik > max_lik){
				maximizers[0]=intercept; 
				maximizers[1]=trend; 
				max_lik = lik ; 
				printf("***");  
			}
			printf("\n");
		}
	}
	// fill with the maximizer 
	fill_trend(params, obssLen, maximizers[0], maximizers[1]); 
	return(max_lik);
}	
	
double find_optimal_trend_greedy
	(double* obss,
	 int obssLen,
	 double* params, // this is the output 
	 char* washSequence,
	 double int_start, 
	 double trend_start,
	 double int_step,
	 double trend_step){
	
	double max_lik,curr_sigma,lik ;
	double intercept;
	double trend; 
	double maximizers[2] ; 
	double res1,res2,res3,res4; 
	int found_better ; 
	int c;
	c=0; 
	intercept=int_start; 
	trend=trend_start ; 
	fill_trend(params, obssLen, intercept, trend); 
	lik = forward(obss,obssLen,params,washSequence);
	while(true){
		c++; 
		fill_trend(params, obssLen, intercept, trend - trend_step); 
		res3 = forward(obss,obssLen,params,washSequence);
		if(res3 > lik){
			lik =res3 ; 
			found_better=true; 
			trend = trend - trend_step ; 
			continue ; 
		}
		fill_trend(params, obssLen, intercept, trend + trend_step); 
		res4 = forward(obss,obssLen,params,washSequence);
		if(res4 > lik){
			lik =res4 ; 
			found_better=true; 
			trend = trend + trend_step ; 
			continue ; 
		}
		found_better=false; 
		fill_trend(params, obssLen, intercept-int_step, trend); 
		res1 = forward(obss,obssLen,params,washSequence);
		if(res1 > lik){
			lik =res1 ; 
			found_better=true; 
			intercept = intercept - int_step ; 
			continue ; 
		}
		fill_trend(params, obssLen, intercept+int_step, trend); 
		res2 = forward(obss,obssLen,params,washSequence);
		if(res2 > lik){
			lik =res2 ; 
			found_better=true; 
			intercept = intercept + int_step ; 
			continue ; 
		}
//		printf("%f %f %f %f %f\n", lik, res1 ,res2 ,res3, res4);
		if(found_better==false){break;}
	}
		
	fill_trend(params, obssLen, intercept, trend); 
	printf("%f\t%f\n",intercept,trend); //now printing this
	//printf("%d\n",c); previously printing this
	return(lik);
}	

//dupliccated the find_optimal_trend_greedy() function and modified it to get the intercept and trend instead of likely hood value.
struct_it return_optimal_trend_intercept_greedy
	(double* obss,
	 int obssLen,
	 double* params, // this is the output 
	 char* washSequence,
	 double int_start, 
	 double trend_start,
	 double int_step,
	 double trend_step
	)
{
	struct_it pass;
	double max_lik,curr_sigma,lik ;
	double intercept;
	double trend; 
	double maximizers[2] ; 
	double res1,res2,res3,res4; 
	int found_better ; 
	int c;
	c=0; 
	intercept=int_start; 
	trend=trend_start ; 
	fill_trend(params, obssLen, intercept, trend); 
	lik = forward(obss,obssLen,params,washSequence);
	while(true){
		c++; 
		fill_trend(params, obssLen, intercept, trend - trend_step); 
		res3 = forward(obss,obssLen,params,washSequence);
		if(res3 > lik){
			lik =res3 ; 
			found_better=true; 
			trend = trend - trend_step ; 
			continue ; 
		}
		fill_trend(params, obssLen, intercept, trend + trend_step); 
		res4 = forward(obss,obssLen,params,washSequence);
		if(res4 > lik){
			lik =res4 ; 
			found_better=true; 
			trend = trend + trend_step ; 
			continue ; 
		}
		found_better=false; 
		fill_trend(params, obssLen, intercept-int_step, trend); 
		res1 = forward(obss,obssLen,params,washSequence);
		if(res1 > lik){
			lik =res1 ; 
			found_better=true; 
			intercept = intercept - int_step ; 
			continue ; 
		}
		fill_trend(params, obssLen, intercept+int_step, trend); 
		res2 = forward(obss,obssLen,params,washSequence);
		if(res2 > lik){
			lik =res2 ; 
			found_better=true; 
			intercept = intercept + int_step ; 
			continue ; 
		}
//		printf("%f %f %f %f %f\n", lik, res1 ,res2 ,res3, res4);
		if(found_better==false){break;}
	}
		
	//fill_trend(params, obssLen, intercept, trend); 
	pass.intercept=intercept;
	pass.trend=trend;
	//printf("%f\t%f\n",intercept,trend); //now printing this
	//printf("%d\n",c); previously printing this
	return(pass);
}	


double find_optimal_trend_only
	(double* obss,
	 int obssLen,
	 double* params, // this is the output 
	 char* washSequence,
	 double intercept, 
	 double trend_start,
	 double trend_end,
	 int steps){
	
	double max_lik,curr_sigma,lik ;
	double step;
	double trend; 
	double maximizer ; 
	step = (trend_end-trend_start)/steps ; 
	
	maximizer=-1; 
	max_lik=-1; 
	for(trend = trend_start ; trend <= trend_end ; trend += step){
		// fill the params vector 
		fill_trend(params, obssLen, intercept,trend); 
		lik = forward(obss,obssLen,params,washSequence);


//		printf("%f-%f\t%f",intercept,trend,lik);

		if(lik > max_lik){
			maximizer=trend; 
			max_lik = lik ; 
//			printf("***");  
		}
//		printf("\n");
	
	}
	// fill with the maximizer 
	fill_trend(params, obssLen, intercept, maximizer); 
	return(max_lik);
}	

	
	//----------
//
// viterbi--
//	Infer a flow sequence from a vector of flowgram observations.
//
//----------
//
// Arguments:
//	double*	obss:			A vector of flowgram observations.
//	int		obssLen:		the number of observations.
//	char*	washSequence:	The periodic sequence of washes corresponding to
//							.. positions in the flow sequence.
//	int*	seqLen:			Place to return the length of the inferred flow
//							.. sequence.
//
// Returns:
//	A pointer to a newly allocated vector containing the inferred flow
//	sequence.  Failures cause program termination.  The caller is responsible
//	for (eventually) disposing of this vector, with a call to free().
//
//----------

static int* viterbi
   (double*		obss,
	int			obssLen,
	double* 	params,
	char*		washSequence,
	int*		_seqLen,
	int	dist)
	{
	introw*		bestPrev = NULL;
	realrow*	logGraph = NULL;
	introw*		incGraph = NULL;
	int*		seq      = NULL;
	int			washLen = strlen(washSequence);
	int			i, j, seqLen, start, curr;

	// allocate the viterbi matrixes
	// nota bene: seq is allocated without what would be entry zero, so all
	//            indexes into it must have 1 subtracted from them

	seqLen = obssLen + 1;

	seq = malloc((seqLen-1)*sizeof(int));
	if (seq == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate seq for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	bestPrev = malloc(seqLen*sizeof(introw));		// transitions
	if (bestPrev == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate bestPrev for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	logGraph = malloc(seqLen*sizeof(realrow));		// likelihoods
	if (logGraph == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate logGraph for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	incGraph = malloc(seqLen*sizeof(introw));		// incorporations
	if (incGraph == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate incGraph for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	// initiate the viterbi matrixes

	for (i=0 ; i<seqLen ; i++)
		{
		for (j=0 ; j<NUM_STATES ; j++) bestPrev[i][j] = 0;
		for (j=0 ; j<NUM_STATES ; j++) logGraph[i][j] = negInf;
		for (j=0 ; j<NUM_STATES ; j++) incGraph[i][j] = 0;
		}

	logGraph[0][NUM_STATES-1] = 0;

	// initiate

	for (i=0 ; i<obssLen ; i++)
		one_step (obss[i], washSequence[i%washLen], i,
		          logGraph[i], logGraph[i+1], bestPrev[i+1], incGraph[i+1],params[i],dist);

	// (post_processing has been moved into the viterbi function)

	start = index_of_max (logGraph[seqLen-1]);

	// trace back through viterbi matrices to build sequence

	i = seqLen-1;
	seq[i-1] = incGraph[i][start];
	curr     = bestPrev[i][start];
	for (i=seqLen-2 ; i>=1 ; i--)
		{
		seq[i-1] = incGraph[i][curr];
		curr     = bestPrev[i][curr];
		}

	if (bestPrev != NULL) free (bestPrev);
	if (logGraph != NULL) free (incGraph);
	if (incGraph != NULL) free (logGraph);

	*_seqLen = seqLen-1;
	return seq;
	}

static int* viterbi_prepros
   (double*	obss,
	int	obssLen,
	double*	params,
	char*	washSequence,
	int*	_seqLen,
	prep_matrix prep,
	char* ml_tag,
	int dist)
	{
	introw*		bestPrev = NULL;
	realrow*	logGraph = NULL;
	introw*		incGraph = NULL;
	int*		seq      = NULL;
	int			washLen = strlen(washSequence);
	int			i, j, seqLen, start, curr;

	// allocate the viterbi matrixes
	// nota bene: seq is allocated without what would be entry zero, so all
	//            indexes into it must have 1 subtracted from them

	seqLen = obssLen + 1;

	seq = malloc((seqLen-1)*sizeof(int));
	if (seq == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate seq for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	bestPrev = malloc(seqLen*sizeof(introw));		// transitions
	if (bestPrev == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate bestPrev for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	logGraph = malloc(seqLen*sizeof(realrow));		// likelihoods
	if (logGraph == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate logGraph for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	incGraph = malloc(seqLen*sizeof(introw));		// incorporations
	if (incGraph == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate incGraph for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}

	// initiate the viterbi matrixes

	for (i=0 ; i<seqLen ; i++)
		{
		for (j=0 ; j<NUM_STATES ; j++) bestPrev[i][j] = 0;
		for (j=0 ; j<NUM_STATES ; j++) logGraph[i][j] = negInf;
		for (j=0 ; j<NUM_STATES ; j++) incGraph[i][j] = 0;
		}

	logGraph[0][NUM_STATES-1] = 0;

	// initiate

	for (i=0 ; i<obssLen ; i++)
		one_step_prepros (obss[i], washSequence[i%washLen], i,
		          logGraph[i], logGraph[i+1], bestPrev[i+1], incGraph[i+1],params[i],prep,ml_tag,dist);

	// (post_processing has been moved into the viterbi function)

	start = index_of_max (logGraph[seqLen-1]);

	// trace back through viterbi matrices to build sequence

	i = seqLen-1;
	seq[i-1] = incGraph[i][start];
	curr     = bestPrev[i][start];
	for (i=seqLen-2 ; i>=1 ; i--)
		{
		seq[i-1] = incGraph[i][curr];
		curr     = bestPrev[i][curr];
		}

	if (bestPrev != NULL) free (bestPrev);
	if (logGraph != NULL) free (incGraph);
	if (incGraph != NULL) free (logGraph);

	*_seqLen = seqLen-1;
	return seq;
	}



static double forward
   (double*		obss,
    int			obssLen,
	double* 	params,
	char*		washSequence)
	{
	realrow*	probGraph = NULL;
	int			washLen = strlen(washSequence);
	int			i, j, seqLen, start, curr;

	// allocate the forward matrixes
	// nota bene: seq is allocated without what would be entry zero, so all
	//            indexes into it must have 1 subtracted from them
	seqLen = obssLen + 1;
	probGraph = malloc(seqLen*sizeof(realrow));		// likelihoods
	if (probGraph == NULL)
		{
		fprintf (stderr, "(in viterbi) Failed to allocate logGraph for n=%d\n", seqLen-1);
		exit (EXIT_FAILURE);
		}


	// initiate the forward matrixes

	for (i=0 ; i<seqLen ; i++)
		{
		for (j=0 ; j<NUM_STATES ; j++) probGraph[i][j] = 0;
		}

	probGraph[0][NUM_STATES-1] = 1;

	// initiate

	for (i=0 ; i<obssLen ; i++)
		one_step_forward (obss[i], washSequence[i%washLen], i,
		          probGraph[i], probGraph[i+1],params[i]);


	// get the sum of the last layer - this is the likelihood
	double lik = 0 ;
	for(i = 0 ; i < NUM_STATES ; i++){
		lik += probGraph[seqLen-1][i] ; 
	}

	if (probGraph != NULL) free (probGraph);

	return lik;
	}
	
//----------
// $$$ verify whether the statement about "uniform distribution of bases" is
//     .. still correct, given that we apparently make use of GC probabilty
//
// create_transition_probs--
//	Create a table TRANS_PROBS with the transition probablities between states
//	based on the assumption of uniform distribution of bases.  For example, if
//	the state is AGC and the wash is C, we switch to AG with probability 1/3
//	(no incorporation) and to AGT with probability 2/3.
//
//----------
//
// Arguments:
//	double	gcContent:	Expected level of GC content in the sequenced sample.
//
// Returns:
//	Nothing, but fills the global array TRANS_PROBS.
//
//----------

static void create_transition_probs
   (double	gcContent)
	{
	double	atContent = 1 - gcContent;
	double	probs[4], newProbs[4], sumNewProbs;
	sstate	oldState, s;
	int		oldStateNum, nextStateNum, wash;
	int		i, j, k;

	//initialize TRANS_PROBS to 0

	for (i=0 ; i<NUM_STATES ; i++)
		for (j=0 ; j<NUM_STATES ; j++)
			for (k=0 ; k<4 ; k++)
				TRANS_PROBS[i][j][k] = 0;

	probs[nt2num('A')] = probs[nt2num('T')] = atContent/2;
	probs[nt2num('C')] = probs[nt2num('G')] = gcContent/2;

	for (oldStateNum=0 ; oldStateNum<NUM_STATES ; oldStateNum++)
		{
		for (wash=0 ; wash<4 ; wash++)
			{
			sstate_from_num (oldState, oldStateNum);
			if (!sstate_get(oldState,wash))
				TRANS_PROBS[oldStateNum][oldStateNum][wash] = 1;
			else
				{
				sumNewProbs = 0;
				for (k=0 ; k<4 ; k++)
					{
					newProbs[k] = 0;
					if (sstate_get(oldState,k))
						{
						newProbs[k] = probs[k]; 
						sumNewProbs += newProbs[k];
						}
					}
				for (k=0 ; k<4 ; k++) newProbs[k] /= sumNewProbs;

				sstate_from_num (s, oldStateNum);
				sstate_skip_num (s, wash);
				nextStateNum = sstate_to_num (s);
				TRANS_PROBS[oldStateNum][nextStateNum][wash] = 1 - newProbs[wash];

				sstate_from_num (s, oldStateNum);
				sstate_see_num  (s, wash);
				nextStateNum = sstate_to_num (s);
				TRANS_PROBS[oldStateNum][nextStateNum][wash] = newProbs[wash];
				}
			}
		}

	}

//----------
//
// index_of_max--
//	Find location in state vector that contains a maximum value.
//
//----------
//
// Arguments:
//	realrow	vals:	The state vector.
//
// Returns:
//	The index (from 0..15) of an entry contaning the maximum value.
//
//----------

static int index_of_max
   (realrow	vals)
	{
	double	maxVal;
	int		maxLoc;
	int		i;

	maxVal = vals[0];
	maxLoc = 0;
	for (i=1 ; i<NUM_STATES ; i++)
		{
		if (vals[i] >= maxVal)
			{ maxVal = vals[i];  maxLoc = i; }
		}

	return maxLoc;
	}

//----------
//
// one_step--
//	Fill one row of the Viterbi matrix.
//
//----------
//
// $$$ TODO describe the interface here
// Arguments:
//	double	obs:
//	char	wash:
//	int		washIndex:
//	realrow	oldScores:
//	realrow	newScores:	
//	introw	newBest:		
//	introw	newBestInc:		
//
// $$$ TODO describe the interface here
// Returns:
//	Nothing, but ...
//
//----------


static void one_step_prepros
   (double	obs,
	char	wash,
	int		washIndex,
	const realrow oldScores,
	realrow	newScores,
	introw	newBest,
	introw	newBestInc,
	double param,
	prep_matrix prep,
	char* ml_tag,
	int dist)
	{
	sstate	s;
	prep_matrix p_mat=prep;
	int		prevState;
	int		nextState[2];
	int		whichMax[2];
	double	p[2];
	double	newScore[2];
	int		j;
	int*** whichMax_matrix=p_mat.whichMax_matrix;
	double*** lik_with_incorp_matrix=p_mat.lik_with_incorp_matrix;
	double** lik_without_incorp_matrix=p_mat.lik_without_incorp_matrix;

	for (prevState=0 ; prevState<NUM_STATES ; prevState++)
		{
		// check if this is a dead lead - a state which is pretty much impossible
		if (oldScores[prevState] <= negInf) continue;

		sstate_from_num  (s, prevState);
		sstate_skip_char (s, wash);
		nextState[0] = sstate_to_num (s);

		sstate_from_num (s, prevState);
		sstate_see_char (s, wash);
		nextState[1] = sstate_to_num (s);

		whichMax[0] = whichMax[1] = 0;
		if(ml_tag=="YES"){
		    //if we directly convert the value we are missing information 1.16 is sometimes coverted to 115 instead of 116	
		   double tmp_std=ceil(obs*1000);
		   int  pos= (int)tmp_std;	
 		   if(lik_without_incorp_matrix[pos/10][washIndex] ){
		   // printf(".");
		    p[0]=lik_without_incorp_matrix[pos/10][washIndex];
		   
 		   }else{
		    p[0] = obs_given_no_incorp (obs,param,dist);
  		   }
		
		  if(lik_with_incorp_matrix[pos/10][nt2num(wash)][washIndex]){
		  //printf("+"); 
		    p[1]=lik_with_incorp_matrix[pos/10][nt2num(wash)][washIndex];
		    //whichMax[1]=whichMax_matrix[(int)(obs*100)][nt2num(wash)][washIndex];
			 int wm=whichMax_matrix[pos/10][nt2num(wash)][washIndex];
			 double tp=obs_given_incorp    (obs, wash, &whichMax[1],param,dist);
			if(p[1]!=tp){
			//printf("%f \t %f\t %f\t %c\t%d\t%d\t%d\t %d\t%d\n",p[1],tp,obs,wash,washIndex,wm,whichMax[1],pos,pos/10);
			}
		  }else{
		    p[1] = obs_given_incorp    (obs, wash, &whichMax[1],param,dist);
		  }
		}else if(ml_tag=="NO"){
			p[0] = obs_given_no_incorp (obs,param,dist);
			p[1] = obs_given_incorp    (obs, wash, &whichMax[1],param,dist);

		}
		p[0] *= TRANS_PROBS[prevState][nextState[0]][nt2num(wash)];
		p[1] *= TRANS_PROBS[prevState][nextState[1]][nt2num(wash)];
		if (p[0] + p[1] == 0) continue;

		for (j=0 ; j<2 ; j++)
			{
			if (nextState[j] == -1) continue;
			newScore[j] = oldScores[prevState] + log(p[j]);
			if (newScore[j] > newScores[nextState[j]])
				{
				newScores [nextState[j]] = newScore[j];
				newBest   [nextState[j]] = prevState;
				newBestInc[nextState[j]] = whichMax[j];
				}
			}
		}
	}


static void one_step
   (double	obs,
	char	wash,
	int		washIndex,
	const realrow oldScores,
	realrow	newScores,
	introw	newBest,
	introw	newBestInc,
	double param,
	int dist)
	{
	sstate	s;
	int		prevState;
	int		nextState[2];
	int		whichMax[2];
	double	p[2];
	double	newScore[2];
	int		j;

	for (prevState=0 ; prevState<NUM_STATES ; prevState++)
		{
		// check if this is a dead lead - a state which is pretty much impossible
		if (oldScores[prevState] <= negInf) continue;

		sstate_from_num  (s, prevState);
		sstate_skip_char (s, wash);
		nextState[0] = sstate_to_num (s);

		sstate_from_num (s, prevState);
		sstate_see_char (s, wash);
		nextState[1] = sstate_to_num (s);

		whichMax[0] = whichMax[1] = 0;
		p[0] = obs_given_no_incorp (obs,param,dist);
		p[1] = obs_given_incorp    (obs, wash, &whichMax[1],param,dist);
		p[0] *= TRANS_PROBS[prevState][nextState[0]][nt2num(wash)];
		p[1] *= TRANS_PROBS[prevState][nextState[1]][nt2num(wash)];
		if (p[0] + p[1] == 0) continue;

		for (j=0 ; j<2 ; j++)
			{
			if (nextState[j] == -1) continue;
			newScore[j] = oldScores[prevState] + log(p[j]);
			if (newScore[j] > newScores[nextState[j]])
				{
				newScores [nextState[j]] = newScore[j];
				newBest   [nextState[j]] = prevState;
				newBestInc[nextState[j]] = whichMax[j];
				}
			}
		}
	}

static void one_step_forward
   (double	obs,
	char	wash,
	int		washIndex,
	const realrow oldProbs,
	realrow	newProbs,
	double param)
	{
	sstate	s;
	int		prevState;
	int		nextState[2];
	double	p[2];
	int		i,j;

	for (prevState=0 ; prevState<NUM_STATES ; prevState++)
		{
		// check if this is a dead lead - a state which is pretty much impossible
		if (oldProbs[prevState] == 0) continue;

		sstate_from_num  (s, prevState);
		sstate_skip_char (s, wash);
		nextState[0] = sstate_to_num (s);

		sstate_from_num (s, prevState);
		sstate_see_char (s, wash);
		nextState[1] = sstate_to_num (s);

		p[0] = obs_given_no_incorp (obs,param,dist);
		p[1] = obs_given_incorp_forward    (obs, wash, washIndex,param,dist);
		p[0] *= TRANS_PROBS[prevState][nextState[0]][nt2num(wash)];
		p[1] *= TRANS_PROBS[prevState][nextState[1]][nt2num(wash)];
		if (p[0] + p[1] == 0) continue;

		for (j=0 ; j<2 ; j++)
			{
			if (nextState[j] == -1) continue;
			newProbs[nextState[j]] += oldProbs[prevState] * p[j];
		}
	}
	for(i=0; i < NUM_STATES ; i++){
		newProbs[i] *= UNDERFLOW_CONST ; 
	}

}

	// probability of observing data given that there was no incorporation;  norm
// is the noise

static double obs_given_no_incorp
   (double	obs,
	double param,
	int dist)
	{
	if(dist==NORMAL){
		return dnorm (obs, 0, param);
	}else{ // double exponential
		return ddexp(obs,0,param); 
	}
}


// use have incoporation
// you don't know how many were incorporated, but for viterbi you want to know
// .. this
// keeping two things, maximum (best value of incorporation) and sum of
// .. likelihoods

static double obs_given_incorp
   (double	obs,
	char	wash,
	int*	whichMax,
	double param,
	int dist)
	{
	double	M;
	double	contentProb;
	double	maxLik;
	double	prob, mean, lik;
	int		i;

	
	M = floor(obs + 1.0) + 2;
	if (M < 4.0) M = 4.0;

	if ((wash == 'A') || (wash == 'T')) contentProb = (1 - gcContent) / 2;
	                               else contentProb = (gcContent) / 2;

	maxLik = -1;
	*whichMax = -1;
	for (i=0 ; i<M ; i++)
		{
		prob = dgeom (i+1, 1 - contentProb); // note that we are conditioning on the fact that we have an incorporation
		mean = i + 1;
		if(dist==NORMAL){
			lik  = prob * dnorm (obs, mean, param);
		}else{ // double exponential
			lik  = prob * ddexp (obs, mean, param);
		}
		
		if (lik >= maxLik)
			{
			maxLik = lik;
			*whichMax = i+1; // $$$ TODO not sure about the indices here...
			}
		}

	return maxLik;
	}


static double obs_given_incorp_forward
   (double	obs,
	char	wash,
	int		washIndex,
	double param,
	int dist)
	{

	double	M;
	double	contentProb;
	double	sumLik=0;
	double	prob, mean, lik;
	int		i;

	M = floor(obs + 1.0) + 2;
	if (M < 4.0) M = 4.0;

	if ((wash == 'A') || (wash == 'T')) contentProb = (1 - gcContent) / 2;
	                               else contentProb = (gcContent) / 2;

	for (i=0 ; i<M ; i++)
		{
		prob = dgeom (i+1, 1 - contentProb); // note that we are conditioning on the fact that we have an incorporation
		mean = i + 1;
		if(dist==NORMAL){
			sumLik  += prob * dnorm (obs, mean, param);
		}else{ // double exponential
			sumLik  += prob * ddexp (obs, mean, param);
		}
	}
	return sumLik;
	}

static double dgeom
   (int		x,
	double	prob)
	{
	// $$$ TODO check if this is correct
	return pow(1 - prob, x-1) * prob;
	}


static double dnorm
   (double	x,
	double	mean,
	double	sd)
	{
	// $$$ TODO check if this is correct, just took from internet
	return (1 / (sd*sqrt(2*M_PI))) * exp(-0.5 * pow((x-mean)/sd,2.0));
	}

	// b is the parameter of the double exponential. 
	// note that using this parameterization the variance
	// is 2b^2
static double ddexp
	(double x,
	 double mean,
	 double sd)
	{
	double b=sqrt(2)*sd;
//	double b=sd/sqrt(2); <- this is the truth, but since it is only a constant difference we do not care
	return ((1/(2*b)) * exp(-(fabs(x-mean)/b)));
	}
//----------
//
// operations on single states (sstate)--
//
// An sstate is a (zero-terminated) four-character string, e.g. "AC-T".  A
// state may also be "invalid", in which case it is an empty string.  There's
// a one-to-one correspondence between these strings and the numbers 0 thru 15.
//
//----------

// initiate state from number 0..15

static void sstate_from_num
   (sstate	str,
	int		statenum)
	{
	strcpy (str, "----");
	if ( statenum      % 2 == 1) str[3] = 'T';
	if ((statenum / 2) % 2 == 1) str[2] = 'G';
	if ((statenum / 4) % 2 == 1) str[1] = 'C';
	if ((statenum / 8) % 2 == 1) str[0] = 'A';
	}

// convert state to number 0..15 (or -1 if the state is invalid)

static int sstate_to_num
   (sstate	str)
	{
	int		statenum = 0;

	if (strlen(str) < 4) return -1;
	if (str[0] != '-') statenum += 8;
	if (str[1] != '-') statenum += 4;
	if (str[2] != '-') statenum += 2;
	if (str[3] != '-') statenum += 1;

	return statenum;
	}

// does state have room to accept wash?

static int sstate_get
   (sstate	str,
	int		wash)
	{
	return (str[wash] != '-');
	}

// change the state like you're incorporating a wash

static void sstate_see_char (sstate str, char c)
	{ sstate_see_num (str, nt2num(c)); }

static void sstate_see_num
   (sstate	str,
	int		x)
	{
	if      (strlen(str) < 4) return;		// (state is already invalid)
	else if (str[x] == '-')   str[0] = 0;	// (make the state invalid)
	else                      { strcpy (str, "ACGT");  str[x] = '-'; }
	}

// change the state like you're not incorporating a wash

static void sstate_skip_char (sstate str, char c)
	{ sstate_skip_num (str, nt2num(c)); }

static void sstate_skip_num (sstate str, int x)
	{ str[x] = '-'; }

// encode a nucleotide (A, C, G, or T) to a number in the range 0..3.

static int nt2num
   (char	c)
   {
	if (c == 'A') return 0;
	if (c == 'C') return 1;
	if (c == 'G') return 2;
	if (c == 'T') return 3;
	if (c == 'a') return 0;
	if (c == 'c') return 1;
	if (c == 'g') return 2;
	if (c == 't') return 3;

	if (isprint(c)) fprintf (stderr, "Invalid character %c in nt2num\n", c);
	           else fprintf (stderr, "Invalid character 0x%02X in nt2num\n", c);

	exit (EXIT_FAILURE);
	}

//----------
//
// open_file--
//	Open a file (this is a proxy for fopen).
//
//----------
//
// Arguments:
//	char*	basename:	The name file to open (but see extension).
//	char*	extension:	If non-NULL, this is appended to the basename to
//							.. for the filename.
//	char*	mode:		File access mode, e.g. "rt", "rb", "wt", or "wb".
//
// Returns:
//	A pointer to the opened file.  Failures cause program termination.  The
//	caller is responsible for disposing of this, with a call to fclose().
//
//----------

static FILE* open_file
   (char*		basename,
	const char*	extension,
	const char*	mode)
	{
	char*	filename = NULL;
	int		nameIsFromHeap = false;
	FILE*	f;
	int		bytesNeeded;

	if (extension == NULL)
		filename = basename;
	else
		{
		bytesNeeded = strlen(basename)
		            + strlen(extension) + 1;
		filename = malloc (bytesNeeded);
		if (filename == NULL)
			{
			fprintf (stderr, "Failed to allocate filename for \"%s%s\"\n",
			         basename, extension);
			exit (EXIT_FAILURE);
			}
		strcpy (filename,                  basename);
		strcpy (filename+strlen(filename), extension);
		}

	f = fopen (filename, mode);
	if (f == NULL)
		{
		fprintf (stderr, "Failed to open \"%s\" for \"%s\"\n",	 filename, mode);
		if (nameIsFromHeap) free (filename);
		exit (EXIT_FAILURE);
		}

	if (nameIsFromHeap) free (filename);

	return f;
	}

//----------
//
// read_line--
//	Read one line from a text file.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to read from.
//	filestate*	fState:		Our internal state associated with that file.  See
//							.. note (1) below.
//	char*		line:		The buffer to receive the line.  See note (2).
//	size_t		lineSize:	The size of the buffer, including room for a
//							terminating zero.
//
// Returns:
//	true if we were successful, false if not (i.e. false at end-of-file).
//	Failures other than end-of-file cause program termination.  See note (2).
//
//----------
//
// Notes:
//
//	(1)	Prior to fState's use by this function, the caller must initialize
//		fState by calling init_read_line.
//
//	(2)	The line is returned with a zero-terminator but without a new line.
//		Checks are made to assure that the line does not exceed the length of
//		buffer.  These checks allow the final line in the file to be without a
//		newline.  If any other line *appears* to have ended without a newline,
//		this indicates the buffer was not large enough, and is considered a
//		failure.  This will not be evident to the caller until the second call
//		for the long line (the first call is returned as successful, but is
//		truncated).
//
//----------

static void init_read_line
   (filestate*	fState,
	char*		filename)
	{
	fState->filename   = filename;
	fState->lineNum    = 0;
	fState->missingEol = false;
	}


static int read_line
   (FILE*		f,
	filestate*	fState,
	char*		line,
	size_t		lineSize)
	{
	int		len;

	while (true)
		{
		if (fgets (line, lineSize, f) == NULL) return false;

		fState->lineNum++;

		// check for lines getting split by fgets (the final line in the file
		// might not have a newline, but no internal lines can be that way)

		if (fState->missingEol)
			{
			fprintf (stderr, "line is too long (%s: line %d)", fState->filename, fState->lineNum-1);
			exit (EXIT_FAILURE);
			}

		len = strlen(line);
		if (len != 0)
			{ fState->missingEol = (line[len-1] != '\n');  break; }
		}

	return true;
	}

//----------
//
// get_observations--
//	Read one line of flowgram observations from a text file.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to read from.
//	filestate*	fState:		The internal state associated with that file, for
//							.. read_line().  See the description of read_line()
//							.. for more details.
//	int*		obssLen:	Place to return the number of observations (the
//							.. length of the observation vector).
//
// Returns:
//	A pointer to a newly allocated vector of observations.  Failures cause
//	program termination.  The caller is responsible for (eventually) disposing
//	of this vector, with a call to free().
//
//----------

static double* get_observations
   (FILE*		f,
	filestate*	fState,
	int*		obssLen)
	{
	static char	line[1000*6+1];	// (room for 1000 floats, up to 5 chars each)
	double*	obss = NULL;
	char*	s, *sBeg, *sEnd;
	int		numFields, ix;

	if (!read_line (f, fState, line, sizeof(line)))
		return NULL;

	// count the number of fields in the line

	s = skip_whitespace (line);
	numFields = 0;
	while (*s != 0)
		{
		s = skip_darkspace  (s);
		s = skip_whitespace (s);
		numFields++;
		}

	// allocate the observation vector

	obss = malloc(numFields*sizeof(double));
	if (obss == NULL)
		{
		fprintf (stderr, "(in get_observations) Failed to allocate obss for %d observations\n",
		         numFields);
		exit (EXIT_FAILURE);
		}

	// parse the fields

	ix = 0;
	s = skip_whitespace (line);
	while (*s != 0)
		{
		sBeg = s;
		sEnd = skip_darkspace (s);
		s = skip_whitespace (sEnd);
		*sEnd = 0;
		obss[ix] = string_to_double (sBeg);
		ix++;
		}

	// success

	*obssLen = numFields;
	return obss;
	}

//----------
//
// zero_trim_observations--
//	Trim zeros from the end of a vector of flowgram observations
//
//----------
//
// Arguments:
//	double*	obss:		The vector to trim.
//	int*	obssLen:	The length of the vector.  If any observations are
//						.. trimmed, this value is changed.
//
// Returns:
//	(nothing)
//
//----------

static void zero_trim_observations
   (double*	obss,
	int*	obssLen)
	{
	int		origLen = (*obssLen);
	int		endIx   = origLen - 1;

	while ((endIx >= 0) && (obss[endIx] == 0)) endIx--;

	if (endIx+3 < origLen) *obssLen = endIx+3;
	                  else *obssLen = origLen;
	}

//----------
//
// write_flow_as_nts--
//	Convert a sequence of flow values to nucleotides and write them to a file.
//
//----------
//
// Arguments:
//	FILE*	f:				The file to write to.
//	int*	seq:			The sequence of flow values to convert.
//	int		seqLen:			The length of the sequence.
//	char*	washSequence:	The periodic sequence of washes corresponding to
//							.. positions in the flow sequence.
//
// Returns:
//	(nothing)
//
//----------

static void write_flow_as_nts
   (FILE*	f,
	int*	seq,
	int		seqLen,
	char*	washSequence)
	{
	int		washLen = strlen(washSequence);
	int		i, j;

	for (i=0 ; i<seqLen ; i++)
		{
		for (j=0 ; j<seq[i] ; j++)
			fprintf (f, "%c", washSequence[i%washLen]);
		}
	fprintf (f, "\n");
	}

//----------
//
// string_to_double--
//	Parse a string for the double floating point value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid number-- failures result in fatality.
//
//----------

static double string_to_double
   (const char*	s)
	{
	double		v;
	char		extra;

	if (sscanf (s, "%lf%c", &v, &extra) != 1)
		{
		fprintf (stderr, "\"%s\" is not a number\n", s);
		exit (EXIT_FAILURE);
		}

	return v;
	}

//----------
//
// skip_whitespace, skip_darkspace--
//	Skip white or dark characters in a zero-terminated string.
//
//----------

static char* skip_whitespace (char* s)
	{ while ((*s != 0) && (isspace (*s))) s++;  return s; }

static char* skip_darkspace (char* s)
	{ while ((*s != 0) && (!isspace (*s))) s++;  return s; }

void print_time()
{
  char buff[100];
  time_t now = time (0);
  strftime (buff, 100, "%Y-%m-%d %H:%M:%S", localtime (&now));
  printf ("Time: %s\n", buff);
  
}
