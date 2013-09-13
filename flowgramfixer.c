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

typedef char sstate[5];		// a string representation of the form "AC-T" or ""

#define NUM_STATES 16
typedef int    introw[NUM_STATES];
typedef double realrow[NUM_STATES];

double TRANS_PROBS[NUM_STATES][NUM_STATES][4];

// state for reading text files line-by-line; see read_line()

typedef struct filestate {
	char*	filename;
	int		lineNum;
	int		missingEol;
	} filestate;

// other globals

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

static void    usage                   (char* programName);


void 		   fill_const			   (double* params, int paramsLen, double fill);

void 		   fill_trend			   (double* params, int paramsLen, double intercept, double trend);
void 		   find_optimal_const_sigma (double* obss,int obssLen,double* params, char* washSequence, double start,	 double end, int steps);
double 		   find_optimal_trend_sigma (double* obss,int obssLen,double* params, char* washSequence, double int_start, double int_end,	 double trend_start, double trend_end, int int_steps, int trend_steps);
double 		   find_optimal_trend_only (double* obss,int obssLen,double* params, char* washSequence, double intercept, double trend_start, double trend_end, int trend_steps);
double 		   find_optimal_trend_greedy (double* obss,int obssLen,double* params, char* washSequence, double int_start, double trend_start, double int_step, double trend_step);

static int*    viterbi                 (double* obss, int obssLen, double* params,
                                        char* washSequence, int* seqLen);
static double  forward                 (double* obss, int obssLen, double* params,
                                        char* washSequence);
static void    create_transition_probs (double gcContent);
static int     index_of_max            (realrow vals);
static void    one_step                (double obs, char wash, int washIndex,
                                        const realrow oldScores,
                                        realrow newScores, introw newBest,
                                        introw newBestInc, double param);
static void    one_step_forward        (double obs, char wash, int washIndex,
                                        const realrow oldProbs,
                                        realrow newProbs, double params);

static double  obs_given_no_incorp     (double obs, char wash, int washIndex,double param,int dist);
static double  obs_given_incorp        (double obs, char wash, int washIndex,
                                        int* whichMax, double param,int dist);
static double  obs_given_incorp_forward(double obs, char wash, int washIndex,double param, int dist);
static double  dnorm                   (double x, double mean, double sd);
static double  ddexp                   (double x, double mean, double sd);
static double  dgeom                   (int x, double prob);
static void    sstate_from_num         (sstate str, int statenum);
static int     sstate_to_num           (sstate str);
static int     sstate_get              (sstate str, int wash);
static void    sstate_see_char         (sstate str, char c);
static void    sstate_see_num          (sstate str, int x);
static void    sstate_skip_char        (sstate str, char c);
static void    sstate_skip_num         (sstate str, int x);
static int     nt2num                  (char c);
static FILE*   open_file               (char* basename, const char* extension,
                                        const char* mode);
static void    init_read_line          (filestate* fState, char* filename);
static int     read_line               (FILE* f, filestate* fState, char* line,
                                        size_t lineSize);
static double* get_observations        (FILE* f, filestate* fState, int* obssLen);
static void    zero_trim_observations  (double* obss, int* obssLen);
static void    write_flow_as_nts       (FILE* f, int* seq, int seqLen,
                                        char* washSequence);
static double  string_to_double        (const char* s);
static char*   skip_whitespace         (char* s);
static char*   skip_darkspace          (char* s);

//----------
//
// main program--
//
//----------

static void usage
   (char* programName)
	{
	fprintf (stderr, "usage: %s incorporation_file output_file dist param1 param2\n", programName);
	fprintf (stderr, "\n");
	fprintf (stderr, "incorporation_file is space delimited file of incorporation sequences\n");
	fprintf (stderr, "output_file is the prefix for all output files\n");
	fprintf (stderr, "dist is the distribution to be used in the likelihood (exp/normal)\n");
	fprintf (stderr, "param1 can be greedy/trend (greedy is greedy search, trend is grid search\n");
	fprintf (stderr, "param2 is an additional parameter for the greedy algorithm (step size - default is 1, larger step size runs faster) \n");
	fprintf (stderr, "alternatively, if param1 and param2 are numbers, they are taken to be the  intercept and trend for the noise model of all the reads\n");
	
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

int main
   (int			argc,
	char**		argv)
	{
	char*		washSequence;
	char*		allObssFilename;
	char*		outputTemplate;
	FILE*		allObssF, *likF, *seqF;
	filestate	allObssState;
	int			readNumber;
	int			i,j;
	double step_size;  
	
	if (argc != 5 && argc != 6) usage (argv[0]);

	allObssFilename = argv[1];
	outputTemplate  = argv[2];
	// parse dist arg 
	if(strcmp(argv[3],"normal")==0){
		dist=NORMAL; 
	}else if(strcmp(argv[3],"exp")==0){
		dist=EXP; 
	}else{
		usage(argv[0]);
	}

	// parse param arg, which gives us the mode
	if(strcmp(argv[4],"intercept")==0){
		mode=INTERCEPT; 
	}else if(strcmp(argv[4],"trend")==0){
		mode=TREND; 
	}else if(strcmp(argv[4],"trendOnly")==0){
		mode=TREND_ONLY;
		sd_start=string_to_double (argv[5]);
	}else if(strcmp(argv[4],"greedy")==0){
		mode=GREEDY ; 
		if(argc==6){
			step_size = string_to_double(argv[5]);
		}else{
			step_size = 1; 
		}
	}else{
		mode=CONST; 
		sd_start  = string_to_double (argv[4]);
		if(argc==6){
			sd_trend  = string_to_double (argv[5]);
		}else{
			sd_trend = 0 ; 
		}
	}
	
	printf("Using the %s distribution (%d) with mode %s %s(%d)\n",argv[3],dist,argv[4],argv[5],mode);

	gcContent = 0.5;
	washSequence = defaultWashSequence;

	// read the additional two arguments 
	// build TRANS_PROBS

	create_transition_probs (gcContent);

	// run viterbi on all strings and dump to files

	allObssF = open_file (allObssFilename, NULL, "rt");
	init_read_line (&allObssState, allObssFilename);
	likF    = open_file (outputTemplate, ".lik", "wt");
	seqF     = open_file (outputTemplate, ".seq",  "wt");

	readNumber = 0;
	while (true)
		{
		double* obss;
		double*    params;
		int*    seq;
		int		obssLen;
		int		seqLen;
		double lik ;
		lik=0;
		
		obss = get_observations (allObssF, &allObssState, &obssLen);
		if (obss == NULL) break;
		readNumber++;
		debugObssLength_1;
		zero_trim_observations (obss, &obssLen);
		debugObssLength_2;

		params = malloc(obssLen * sizeof(double));
		if(mode==CONST){
			fill_trend(params,obssLen,sd_start,sd_trend) ; 
			seq = viterbi (obss, obssLen, params, washSequence, &seqLen);
		}else if(mode==INTERCEPT){
			find_optimal_const_sigma(obss, obssLen, params, washSequence,0,0.2,10);
			seq = viterbi (obss, obssLen, params, washSequence, &seqLen);
		}else if(mode==TREND_ONLY){
			lik=find_optimal_trend_only(obss, obssLen, params, washSequence,sd_start,0,0.0005,10);
			seq = viterbi (obss, obssLen, params, washSequence, &seqLen);
		}else if(mode==GREEDY){
			lik=find_optimal_trend_greedy(obss, obssLen, params, washSequence,0.03,0.0002,0.001*step_size,0.00001*step_size);
			seq = viterbi (obss, obssLen, params, washSequence, &seqLen);
		}else{ // mode is trend 
			lik=find_optimal_trend_sigma(obss, obssLen, params, washSequence,0,0.1,0,0.0005,10,20);
			seq = viterbi (obss, obssLen, params, washSequence, &seqLen);
		}
		fprintf(likF,"%f\t%f\t%f\n",params[0],params[1]-params[0],lik);
			debugSeqLength_1

		// write flow to file

/*		for (j=0 ; j<seqLen-1 ; j++)
			fprintf (flowF, "%d,", seq[j]);
		fprintf (flowF, "%d\n", seq[seqLen-1]);
*/
		// write string to file

		write_flow_as_nts (seqF, seq, seqLen, washSequence);
		free (params);
		free (obss);
		free (seq);
		fflush(likF);
		fflush(seqF);
		}

	
	fclose (allObssF);
	fclose (likF);
	fclose (seqF);

	return EXIT_SUCCESS;
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
	printf("%d\n",c);
	return(lik);
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
	int*		_seqLen)
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
		          logGraph[i], logGraph[i+1], bestPrev[i+1], incGraph[i+1],params[i]);

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


static void one_step
   (double	obs,
	char	wash,
	int		washIndex,
	const realrow oldScores,
	realrow	newScores,
	introw	newBest,
	introw	newBestInc,
	double param)
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
		p[0] = obs_given_no_incorp (obs, wash, washIndex,param,dist);
		p[1] = obs_given_incorp    (obs, wash, washIndex, &whichMax[1],param,dist);
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

		p[0] = obs_given_no_incorp (obs, wash, washIndex,param,dist);
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
	char	wash,
	int		washIndex,
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
	int		washIndex,
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
		fprintf (stderr, "Failed to open \"%s\" for \"%s%s\"\n",
				 filename, mode);
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
