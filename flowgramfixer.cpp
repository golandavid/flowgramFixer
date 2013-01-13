#include<cassert>
#include<cmath>
#include<string>
#include<stdio.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<vector>
#include<map>
#include<set>
#include<list>
#include<queue>
#include<cstdarg>
#include<algorithm>
#include<limits>
//#include<random>

using namespace std;
#include "defs.h"



double Inf = std::numeric_limits<double>::max(); 
double TRANS_PROBS[16][16][4];
double gc_content;
double global_sd_start; 
double global_sd_end; 
double global_sd; 
int UNDERFLOW_CONST=2;


int TRANS_INC[4][16]; 
int TRANS_NO_INC[4][16]; 

#define PRIOR_TABLE_SIZE 20
double BAYES_PRIOR[4][PRIOR_TABLE_SIZE];

//class for a single state (SState)
//
class SState {
	public:
		string str; //a string representation of the form "AC-T"
		bool invalid; 

		//initiate null state
		SState () : invalid(false) { };

		//initiat state from string
		SState (string _str) : str(_str), invalid(false) { }

		//initiate state from number 0..15
		SState (int statenum) : invalid(false) {
			str = "----";
			if (statenum % 2 == 1) str[3] = 'T';
			if ((statenum / 2) % 2 == 1) str[2] = 'G';
			if ((statenum / 4) % 2 == 1) str[1] = 'C';
			if ((statenum / 8) % 2 == 1) str[0] = 'A';
		}

		//get state number
		int num () {
			if (invalid || (str.length() < 4)) return -1;
			int statenum = 0;
			if (str[0] != '-') statenum += 8;
			if (str[1] != '-') statenum += 4;
			if (str[2] != '-') statenum += 2;
			if (str[3] != '-') statenum += 1;
			return statenum;
		}

		//does state have room to accept wash?
		bool get (int wash) { return str[wash] != '-'; }

		//change the state like you're incorporating a wash
		SState * see (char c) {
			return see(nt2num(c)); 
		}
		SState * see (int x) {
			if (invalid || str[x] == '-') {
				invalid = true;
			} else {
				str = "ACGT";
				str[x] = '-';
			}
			return this;
		}

		//change the state like you're not incorporating a wash
		SState * skip (char c) {
			return skip(nt2num(c));
		}

		SState * skip (int x) {
			str[x] = '-';
			return this;
		}
};



//Find location in vector that contains a maximum value 
int index_of_max(vector<double> & vals) {
	assert(vals.size() > 0);
	double max_val = vals.at(0);
	int max_loc = 0;

	for (int i = 1; i < vals.size(); i++) {
		if (vals.at(i) >= max_val) {
			max_val = vals.at(i);
			max_loc = i;
		}
	}
	return max_loc;
}



double dgeom(int x, double prob) {
	//TODO check if this is correct
	return pow (1 - prob, x-1) * prob;
}

double dnorm(double x, double mean, double sd) {
	//TODO check if this is correct, just took from internet
	return ( 1 / ( sd * sqrt(2*M_PI) ) ) * exp( -0.5 * pow( (x-mean)/sd, 2.0 ) );
}

// prob of observing data given that there was no incorporation
// norm is the noise
double obs_given_no_incorp (double obs, char wash, int wash_index, double sd) {
	return dnorm(obs, 0, sd); 
}

//use have incoporation
//you don't know how many were incorporated, but for viterbi you want to know this
//keeping two things, maximum (best value of incorporation) and sum of likelihoods

double obs_given_incorp(double obs, char wash, int wash_index, int & which_max, double sd) {
	double content_prob;
	if ((wash == 'a') || (wash == 't')) {
		content_prob = (1 - gc_content) / 2;
	} else {
		content_prob = (gc_content) / 2;
	}
	// make sure we do not go over the prior table size 
	double M = max(round(obs + 0.5) + 2, 4.0);
	if(M>PRIOR_TABLE_SIZE){M=PRIOR_TABLE_SIZE;}

	double max_lik = -1;
	which_max = -1;
	for (int i = 0; i < M; i++) {
//		double prob = dgeom(i+1, 1 - content_prob); // note that we are conditioning on the fact that we have an incorporation
		double prob = BAYES_PRIOR[nt2num(wash)][i];
		double mean = i + 1;
		double lik = prob * dnorm(obs, mean, sd);
		if (lik >= max_lik) {
			max_lik = lik;
			which_max = i+1; //TODO not sure about the indices here...
		}
	}
	return max_lik;
}

double obs_given_incorp_forward(double obs, char wash, int wash_index, double sd) {
	double content_prob;
	if ((wash == 'a') || (wash == 't')) {
		content_prob = (1 - gc_content) / 2;
	} else {
		content_prob = (gc_content) / 2;
	}
	double M=4;
	if(obs>2){ M=2+round(obs+0.51);}
//	double M = max(round(obs + 0.5) + 3, 4.0);

	double sum_lik = 0;
	for (int i = 0; i < M; i++) {
		double prob = dgeom(i+1, 1 - content_prob); // note that we are conditioning on the fact that we have an incorporation
		double mean = i + 1;
		sum_lik += prob * dnorm(obs, mean, sd);
	}
	return sum_lik;
}


// create a table TRANS_PROBS with the transition probablities between states
// based on the assumption of uniform distribution of bases.
// for example, if the state is AGC and the wash is C, we switch to
// AG with probability 1/3 (no incorporation) and to AGT with probability 2/3.

void create_transition_probs(double gc_content = 0.5) {
	//initialize TRANS_PROBS to 0
	for (int i = 0; i < 16; i++) 
		for (int j = 0; j < 16; j++) 
			for (int k = 0; k < 4; k++) 
				TRANS_PROBS[i][j][k] = 0;

	double at_content = 1 - gc_content;
	double probs[4];
	probs[0] = at_content/2;
	probs[1] = gc_content/2;
	probs[2] = gc_content/2;
	probs[3] = at_content/2;

	for (int old_state_idx = 0; old_state_idx < 16; old_state_idx++) {
		for (int wash = 0; wash < 4; wash++) {
			SState old_state(old_state_idx);
			if (old_state.get(wash)) {
				double new_probs[4];
				double sum_new_probs = 0;
				for (int i = 0; i < 4; i++) {
					new_probs[i] = probs[i] * int(old_state.get(i));
					sum_new_probs += new_probs[i];
				}
				for (int i = 0; i < 4; i++)  new_probs[i] /= sum_new_probs;

				TRANS_PROBS[old_state_idx][SState(old_state_idx).skip(wash)->num()][wash] = 1 - new_probs[wash];
				TRANS_PROBS[old_state_idx][SState(old_state_idx).see (wash)->num()][wash] = new_probs[wash];
			} else {
				TRANS_PROBS[old_state_idx][old_state_idx][wash] = 1;
			}
		}
	}
	return;
}

void create_transitions(){
	for (int i = 0; i < 16; i++){ 
		for (int wash = 0; wash < 4; wash++){ 
			TRANS_INC[wash][i] = SState(i).see(wash)->num() ; 
			TRANS_NO_INC[wash][i] = SState(i).skip(wash)->num() ; 
		}
	}
	return;
}

void init_priors(){
	gc_content=0.5; 
	for(int nuc=0 ; nuc<4 ;nuc++){
		for(int mer=0; mer < 20 ; mer++){
			BAYES_PRIOR[nuc][mer] = dgeom(mer+1, 1 - gc_content/2);
		}
	}
	return;
}

void one_step(double obs, char wash, int wash_index, const vector<double> & old_scores, 
		vector<double> & new_scores, vector<int> & new_best, vector<int> & new_best_inc, double sd) {
/*
	if(wash_index <225 && wash_index >220){
		for(int i=0 ; i < 16; i++){
			if(old_scores[i] == -1 * Inf){
				cout << 'N' << "'" ;
			}else{
				cout << old_scores[i] << ',';
			}
		}
	cout << endl;
	}
*/

		for (int prev_state = 0; prev_state < 16; prev_state++) {
		// check if this is a dead lead - a state which is pretty much impossible
		if (old_scores.at(prev_state) <= -1 * Inf)  continue;

		int next_state[2];
		next_state[0] = SState(prev_state).skip(wash)->num();
		next_state[1] = SState(prev_state).see(wash)->num();

		int which_max[2] = {0,0};
		double p[2];
		p[0] = obs_given_no_incorp(obs, wash, wash_index, sd);
		p[1] = obs_given_incorp(obs, wash, wash_index, which_max[1], sd);
		p[0] *= TRANS_PROBS[prev_state][next_state[0]][nt2num(wash)]; 
		p[1] *= TRANS_PROBS[prev_state][next_state[1]][nt2num(wash)]; 
		if (p[0] + p[1] == 0) continue;

		double new_score[2];
		for (int j = 0; j < 2; j++) {
			if (next_state[j] == -1) continue;
			new_score[j] = old_scores.at(prev_state) + log(p[j]);
			if (new_score[j] > new_scores.at(next_state[j])) {
				new_scores.at(next_state[j]) = new_score[j];
				new_best.at(next_state[j]) = prev_state;
				new_best_inc.at(next_state[j]) = which_max[j];
			} 
		}
	}
	return;
}

void one_step_forward(double obs, char wash, int wash_index, const vector<double> & old_probs, 
		vector<double> & new_probs, double sd) {
		int top=15;
		if(wash_index==0){top=16;}
		for (int prev_state = 1; prev_state < top; prev_state++) {
		// check if this is a dead lead - a state which is pretty much impossible
		// currently this is left out but should return for faster performance
		if (old_probs.at(prev_state) == 0 )  continue;
	
		int next_state[2];
		next_state[0] = TRANS_NO_INC[nt2num(wash)][prev_state]; 
		next_state[1] = TRANS_INC[nt2num(wash)][prev_state]; 
//		next_state[0] = SState(prev_state).skip(wash)->num();
//		next_state[1] = SState(prev_state).see(wash)->num();

		double p[2];
		p[0] = obs_given_no_incorp(obs, wash, wash_index,sd );
		p[1] = obs_given_incorp_forward(obs, wash, wash_index, sd );
		p[0] *= TRANS_PROBS[prev_state][next_state[0]][nt2num(wash)]; 
		p[1] *= TRANS_PROBS[prev_state][next_state[1]][nt2num(wash)]; 
		if (p[0] + p[1] == 0) continue;

		for (int j = 0; j < 2; j++) {
			if (next_state[j] == -1) continue;
			new_probs[next_state[j]] += old_probs.at(prev_state)  * p[j] ; 
		}
	}
	for(int i=0; i < new_probs.size() ; i++){
		new_probs[i] *= UNDERFLOW_CONST ; 
	}
	return;
}


vector<int> viterbi(vector<double> obss, string washs, vector<double> sds) {

	//initiate the viterbi matrixes
	vector<int> temp_row1(16, 0);
	vector<double> temp_row2(16, -1 * Inf);
	vector< vector <int> > best_prev(obss.size() + 1, temp_row1); //transitions
	vector< vector <double> > log_graph(obss.size() + 1, temp_row2); //likelihoods
	vector< vector <int> > inc_graph(obss.size() + 1, temp_row1);  //incorporations


	log_graph.at(0).at(15) = 0;

	// initiate 
	for (int i = 0; i < obss.size(); i++) {
		one_step(obss.at(i), washs.at(i % washs.length()), i, log_graph.at(i), log_graph.at(i+1), best_prev.at(i+1), inc_graph.at(i+1), sds.at(i));
	}

	//I moved the post_processing into the viterbi function
	int n = obss.size() + 1;
	int start = index_of_max(log_graph.back());

	//trace back through viterbi matrices to build sequence
	vector<int> seq(n);
	seq.back() = inc_graph.back().at(start);
	int curr = best_prev.back().at(start); 
	for (int i = seq.size() - 2; i >= 1; i--) {
		seq.at(i) = inc_graph.at(i).at(curr);
		curr = best_prev.at(i).at(curr);
	}
	
	seq.erase(seq.begin()); 
	return(seq); 
}	

double forward(vector<double> obss, string washs, vector<double> sds) {

	//initiate the forward matrixes
	vector<double> temp_row1(16, 0);
	vector<vector<double> > forward_graph(obss.size() + 1, temp_row1); //likelihoods
	// set the first state
	forward_graph.at(0).at(15) = 1;
	
	// initiate 
	for (int i = 0; i < obss.size(); i++) {
/*		for(int j = 0 ; j < 16 ; j++){
			cout << i << ":" << forward_graph.at(i).at(j) << " " ; 
		}
		cout << endl; */
	one_step_forward(obss.at(i), washs.at(i % washs.length()), i, forward_graph.at(i), forward_graph.at(i+1), sds.at(i));
	}

	// return the sum of the last row
	double lik = 0 ; 
	for( int i = 0 ; i < forward_graph.at(obss.size()).size() ; i++){
		lik += forward_graph.at(obss.size()).at(i) ; 
	}
	return(lik); 
}	

vector<double> find_optimal_const_sigma(vector<double> obss, string washs, double start, double end, int steps){
	
	double max_lik = -1 ; 
	double maximizer = start ; 
	for(double curr_sigma = start ; curr_sigma <= end ; curr_sigma += (end-start)/steps){
		vector<double> sds(obss.size(),curr_sigma);
		double lik = forward(obss,washs,sds); 
		cout << curr_sigma << "\t" << lik ; 
		if(lik > max_lik){
			maximizer=curr_sigma; 
			max_lik = lik ; 
			cout << "***";  
		}
		cout << endl;
	}
	vector<double> ret(obss.size(), maximizer); 
	return(ret);
}	

vector<double> find_optimal_trend_sigma(vector<double> obss, string washs, double start1, double start2, double end1, double end2, int steps){
	

	double max_lik = -1 ; 
	double maximizer[2] = {-1,-1} ; 
	double step1 = (end1-start1)/steps;
	double step2 = (end2-start2)/steps;
	for(double sig_start = start1 ; sig_start <= end1 ; sig_start += step1){
		for(double sig_end = start2 ; sig_end <= end2 ; sig_end += step2){
			double inc = (sig_end-sig_start)/obss.size() ; 
			vector<double> sds(obss.size(),0);
			for(int i = 0 ; i < obss.size() ; i++){
				sds[i] = sig_start + inc * i ; 
			}
			
			double lik = forward(obss,washs,sds); 
//			cout << sig_start << "\t" << sig_end << "\t" << lik ; 
			if(lik > max_lik){
				maximizer[0]=sig_start; 
				maximizer[1]=sig_end; 
				max_lik = lik ; 
//				cout << "***";  
			}
//			cout << endl;
		}
	}
	cout << maximizer[0] << "\t" << maximizer[1] << "\t" << max_lik << endl;
	vector<double> sds_final(obss.size(),0);
	double inc = (maximizer[1] - maximizer[0])/obss.size() ; 
	for(int i = 0 ; i < obss.size() ; i++){
			sds_final[i] = maximizer[0] + inc * i ; 
	}
	
	return(sds_final);
}	


string flow_to_seq(vector<int> seq, string washs){
	//convert string of incorporations to string of nts
	//cout << seq << endl;
	string retval = "";
	for (int i = 0; i < seq.size(); i++) {
		for (int j = 0; j < seq.at(i); j++) {
			retval += washs.at(i % washs.size());
		}
	}
	return retval;
}


void update_priors(string filename, string washs){
	// init counters
	unsigned int gc_counter = 0 ;
	unsigned int gc_total = 0 ; 
	double mer_counters[4][PRIOR_TABLE_SIZE] ; 
	unsigned int mer_total[4] = {0,0,0,0} ;
	for(int i = 0 ; i < 4 ; i++){
		for(int j = 0 ; j < PRIOR_TABLE_SIZE ; j++){
			mer_counters[i][j] = 0 ; 
		}
	}

	// open file, use rounding to call sequences, and compute the various priors.
	ifstream all_obss_file;
	open_file(all_obss_file, filename);


	vector<string> row;
	while (get_row(all_obss_file, row, ' ')) {
		for (int i = 0; i < row.size(); i++){
			int curr=round(atof(row[i].c_str()));
			if(curr==0){continue;}
			char curr_wash = washs.at(i%washs.length()); 
			int curr_nuc = nt2num(curr_wash) ; 
			// update GC content counters
			if(curr_wash=='c' || curr_wash=='g'){
				gc_counter += curr;
			}
			gc_total += curr ; 
			// update homopolymer counters
			mer_counters[curr_nuc][min(curr,PRIOR_TABLE_SIZE)-1]++; 
			mer_total[curr_nuc]++;
		}
	}
	// compute actual priors 
	cout << "gc before:" << gc_content << endl;
	gc_content = double(gc_counter)/double(gc_total);
	cout << "gc after:" << gc_content << endl;
	for(int i = 0 ; i < 4 ; i++){
		for(int j = 0 ; j < PRIOR_TABLE_SIZE ; j++){
			BAYES_PRIOR[i][j] = double(mer_counters[i][j])/double(mer_total[i]) ; 
		}
	}

	return; 
}
void write_priors_to_file(string filename){
	ofstream of ; 
	open_file(of, filename);
	of << "GC content:" << gc_content << endl;;
	of << "prior table" << endl; 
	for(int i = 0 ; i < 4 ; i++){
		for(int j = 0 ; j < PRIOR_TABLE_SIZE ; j++){
			of << BAYES_PRIOR[i][j] << "\t" ;
		}
		of << endl; 
	}
	return; 
}
	

void usage(int argc, char * argv[]) {
	cerr << "Usage: " << argv[0] << " incorporation_file output_file [mode]" << endl;
	cerr << "Program descrption:\nincorporation_file is space delimited file of incorporation sequences" << endl;
	cerr << "output_file is the prefix for all output files \n" << endl;
	cerr << "mode is one of three modes currently supported: \n" << endl;
	cerr << "trend fits a trend and intercept for every read \n" << endl;
	cerr << "intercept fits an intercept for every read \n" << endl;
	cerr << "otherwise you need to specify sigma parameter (recommended value: 0.3) \n" << endl;
	exit(1);
}

void trim_obs(vector<double>& obs){
	int end=obs.size()-1 ; 
	int orig=end+1;  
	while(obs.at(end)==0){
		end-- ; 
	}
	obs.resize(min(end+3,orig));
}


int main(int argc, char * argv[]) {

	if ( (argc<4) ) usage(argc,argv);

	string washs = "tacgtacgtctgagcatcgatcgatgtacagc";
	init_priors(); 

	string output_filename(argv[2]); 
	bool trend,intercept,constant;
	double global_sd; 
	// parse mode argument 
	if(string(argv[3])=="trend"){
		cout << "trend mode" <<endl;
		trend=true;
		intercept=false;
		constant=false; 
	}else if(string(argv[3])=="intercept"){
		cout << "intercept mode" <<endl; 
		trend=constant=false;
		intercept=true; 
	}else if(0!=atof(argv[3])){
		trend = intercept = false;
		constant=true;
		global_sd = atof(argv[3]);
		cout << "fixed at " << global_sd << endl;
	}else{
		usage(argc,argv); 
	}
	if(argc==5 && string(argv[4])=="bayes"){
		cout << "Empirical Bayes" << endl; 
		update_priors(string(argv[1]),washs); 
	}
	else{
		cout << "no bayes" << endl; 
	}
	write_priors_to_file(string(argv[2]) + ".priors");
	
	//build TRANS_PROBS 
	create_transition_probs(gc_content);
	create_transitions(); 
	//run viterbi on all strings and dump to files 
	ofstream lik_file, seq_file ; 
	open_file(lik_file , output_filename + ".lik"); 
	open_file(seq_file , output_filename + ".seq"); 

	string all_obss_filename = argv[1];
	ifstream all_obss_file;
	open_file(all_obss_file, all_obss_filename);
	vector<string> row;
	while (get_row(all_obss_file, row, ' ')) {
		vector<double> obs(row.size());
		for (int i = 0; i < row.size(); i++) obs[i] = atof(row[i].c_str());

		trim_obs(obs); 
//		cout << forward(obs,washs) << endl; 
		vector<double> sds;
		if(trend){
			sds = find_optimal_trend_sigma(obs, washs, 0.05,0.18,0.18,0.35, 5);
			lik_file << sds.at(0) << "\t" << sds.at(obs.size()-1) << endl; 
		}else if(intercept){
			sds = find_optimal_const_sigma(obs, washs, 0.01,0.35, 10);
			lik_file << sds.at(0) << endl; 
		}else{
			sds = vector<double>(obs.size(),global_sd); 
			lik_file << sds.at(0) << endl; 
		}
		vector<int> seq = viterbi (obs, washs, sds);
		// write flow to file 
//		for(int j=0; j < seq.size()-1 ; j++){
//			flow_file << seq.at(j) << "," ; 
//		}
//		flow_file << seq.at(seq.size()-1) << endl; 
		
		string s = flow_to_seq(seq,washs);
		// write string to file 
		seq_file << s << endl;
	}
	return 0;
}

