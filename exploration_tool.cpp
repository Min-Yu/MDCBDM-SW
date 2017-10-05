#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include "memory.h"
#include "exploration_tool.h"
#include "pthread.h"

using namespace std;

struct StageTwoCompactionDataSet {
	short int **StageTwoCompactionData;
};

struct VisitTableSet {
	int **VisitTableData;
	short int *VisitTableCycleMapToSigData;
};
#ifdef MULTITHREADING
void *MSCM(void *);
#else
void MSCM(benchmark_listnode);
#endif

void MSCM_for_decision_mode(benchmark_listnode,int ,short int);

void ResultRecordCreation(benchmark_linkedlist ,int ,int);

void StageOneCompaction(short int ,short int*,short int* ); 

void VisitTableCreation(short int, int, VisitTableSet *);
void RowVisitTableCreation(int,VisitTableSet *);
void ColVisitTableCreation(int,VisitTableSet *);
void DiaVisitTableCreation(int,VisitTableSet *);
void RevDiaVisitTableCreation(int,VisitTableSet *);
void RanVisitTableCreation(int,VisitTableSet *);

void StageTwoXoredCompaction(int, int **, short int**, short int**);
void StageTwoMISRCompaction(int, int **, short int**, short int**);
void StageTwoCompaction(short int ,int ,short int ** ,VisitTableSet *  ,StageTwoCompactionDataSet *);
void Comparison(short int,int,StageTwoCompactionDataSet *,StageTwoCompactionDataSet *,bool **);
void Intersection(short int,int, VisitTableSet *, bool **, bool **);

bool ValueCmpForIntSmallToBig(int const & a, int const & b)
{
    return a < b;
}

short int STAGE_ONE_COMPACTION_DATA_SIZE;
short int STAGE_TWO_COMPACTION_DATA_SIZE;

int main(void) {
	srand(time(NULL));
	string iscas85_circuit_list[] = {"c880", "c1908", "c2670", "c432", "c1355", "c499", "c5315", "c6288", "c7552", "c3540"};
	string boom_circuit_list[] = {"spmv.riscvBOOMTile"};
	string circuit_group_list[] = {"for_iscas85_circuits", "for_boom_tilelevel", "for_boom_corelevel", "for_boom_componentlevel" , "for_boom_all_levels", "for_all_circuits"};
	string selected_benchmark_list[MAXIMUM_CYCLES];
	int i,j;
	int selected_mode, win_size_max, win_size_min, fault_num, ran_round;
	float er_wac_one, er_wac_two, er_wac_odd, er_wac_even, er_bc;
	double error_rate[5];	
	//0: error rate between cycles
	//1: error rate within a cycle (one-bit)
	//2: error rate within a cycle (two-bit) 
	//3: error rate within a cycle (odd-bit) 
	//4: error rate within a cycle (even-bit) 
	int selected_benchmark, num_benchmakrs;
	cout << "Multiple Signature Compaction Method Estimation Tool" << endl;
#ifndef FIXED_INPUTS 
	cout << "0: Evaluation mode\t1:Decision mode" << endl;
	cout << "Enter the mode you prefer to use: " ;
	cin >> selected_mode;
	cout << "----------ISCAS'85 circuits----------" << endl;
	for(i=0;i<ISCAS85_CIRCUITS;i++) {
		cout << i+1 << ") " << iscas85_circuit_list[i] << " \t";
		if((i+1)%4 == 0)
			cout << endl;
	}
	cout << endl;
	cout << "-------------------------------------" << endl;
	cout << "--------RISCV BOOM processor---------" << endl;
	for(i=0;i<BOOM_CIRCUITS;i++) {
		cout << i+1+BOOM_CIRCUITS_INDEX << ") " << boom_circuit_list[i] << "\t";
		if((i+1)%3 == 0)
			cout << endl;
	}
	cout << "-------------------------------------" << endl;
	//cout << "-------------------------------------" << endl;
	//cout << "----------Customer circuits----------" << endl;
	//cout << "-------------------------------------" << endl;
	cout << "-------------------------------------" << endl;
	cout << "---------For circuit groups----------" << endl;
	for(i=0;i<CIRCUIT_GROUPS;i++) {
		cout << i+1+CIRCUIT_GROUPS_INDEX << ") " << circuit_group_list[i] << endl;
	}
	cout << "-------------------------------------" << endl;
	cout << "Enter the benchmark you prefer to use: " ;
	cin >> selected_benchmark;
	cout << endl << "----Enter the window size you prefer to observe: ----" << endl;
	cout << "Enter the minimum window size you prefer to observe: ";
	cin >> win_size_min;
	cout << "Enter the maximum window size you prefer to observe: ";
	cin >> win_size_max;
	cout << "-----------------------------------------------------" << endl;
	cout << endl << "-------------Enter the error rate you prefer to generate errors: -------------" << endl;
	cout << "Enter the error rate within a cycle (one-bit) you prefer to generate errors: ";
	cin >> error_rate[1];
	cout << "Enter the error rate within a cycle (two-bit) you prefer to generate errors: ";
	cin >> error_rate[2];
	cout << "Enter the error rate within a cycle (odd-bit) you prefer to generate errors: ";
	cin >> error_rate[3];
	cout << "Enter the error rate within a cycle (even-bit) you prefer to generate errors: ";
	cin >> error_rate[4];
	cout << "Enter the error rate between cycles you prefer to generate errors: ";
	cin >> error_rate[0];
	cout << "------------------------------------------------------------------------------" << endl;
	cout << "Enter how many times of simulation you would like to run: ";
	cin >> ran_round;
#else
	cout << "0: Evaluation mode\t1:Decision mode" << endl;
	cout << "Enter the mode you prefer to use: 1" << endl;
	selected_mode = 1;
	cout << "----------ISCAS'85 circuits----------" << endl;
	for(i=0;i<ISCAS85_CIRCUITS;i++) {
		cout << i+1 << ") " << iscas85_circuit_list[i] << " \t";
		if((i+1)%4 == 0)
			cout << endl;
	}
	cout << endl;
	cout << "-------------------------------------" << endl;
	cout << "--------RISCV BOOM processor---------" << endl;
	for(i=0;i<BOOM_CIRCUITS;i++) {
		cout << i+1+BOOM_CIRCUITS_INDEX << ") " << boom_circuit_list[i] << "\t";
		if((i+1)%3 == 0)
			cout << endl;
	}
	cout << "-------------------------------------" << endl;
	//cout << "-------------------------------------" << endl;
	//cout << "----------Customer circuits----------" << endl;
	//cout << "-------------------------------------" << endl;
	cout << "-------------------------------------" << endl;
	cout << "---------For circuit groups----------" << endl;
	for(i=0;i<CIRCUIT_GROUPS;i++) {
		cout << i+1+CIRCUIT_GROUPS_INDEX << ") " << circuit_group_list[i] << endl;
	}
	cout << "-------------------------------------" << endl;
	cout << "Enter the benchmark you prefer to use: 12" << endl ;
	selected_benchmark = 12;
	cout << endl << "----Enter the window size you prefer to observe: ----" << endl;
	cout << "Enter the minimum window size you prefer to observe: 64" << endl;
	win_size_min = 64;
	cout << "Enter the maximum window size you prefer to observe: 16384" << endl;
	win_size_max = 16384;
	cout << "-----------------------------------------------------" << endl;
	cout << endl << "-------------Enter the error rate you prefer to generate errors: -------------" << endl;
	cout << "Enter the error rate within a cycle (one-bit) you prefer to generate errors: 25" << endl;
	error_rate[1] = 25;
	cout << "Enter the error rate within a cycle (two-bit) you prefer to generate errors: 25" << endl;
	error_rate[2] = 25;
	cout << "Enter the error rate within a cycle (odd-bit) you prefer to generate errors: 25" << endl;
	error_rate[3] = 25;
	cout << "Enter the error rate within a cycle (even-bit) you prefer to generate errors: 25" << endl;
	error_rate[4] = 25;
	cout << "Enter the error rate between cycles you prefer to generate errors: 0.1" << endl;
	error_rate[0] = 0.1;
	cout << "------------------------------------------------------------------------------" << endl;
	cout << "Enter how many times of simulation you would like to run: 100" << endl;
	ran_round = 100;
#endif
	cout << endl << "Initializing..." << endl;
	//Extraction Benchmarks information
	benchmark_linkedlist list;
	benchmark_listnode pending_benchmark;
	ifstream infile("./signal_trace_data/Extracted_signals_cycles.txt");
	string b_name;
	int b_s, b_c;

	//Ensure range of window size
	int sqrt_win_size_max, sqrt_win_size_min, num_windows;
	sqrt_win_size_min = sqrt(win_size_min);
	sqrt_win_size_max = sqrt(win_size_max);
	num_windows = (log2(sqrt_win_size_max) - log2(sqrt_win_size_min) + 1);

	int *sqrt_win_array = new int[num_windows];
	for(i=0;i<num_windows;i++) {
		sqrt_win_array[i] = sqrt_win_size_min * pow(2,i);
#ifdef DEBUG_EN
		cout << "sqrt_win_array[" << i << "] = " << sqrt_win_array[i] << endl;
#endif	
	}
#ifdef DEBUG_EN
	cout << "window size extraction information" << endl;
	cout << "sqrt_win_size_min = " << sqrt_win_size_min << endl;
	cout << "sqrt_win_size_max = " << sqrt_win_size_max << endl;
	cout << "num_windows = " << num_windows << endl;
#endif

	switch (selected_benchmark) {
		case CIRCUIT_GROUPS_ISCAS85:
			for(i=0;i<ISCAS85_CIRCUITS;i++) {
				selected_benchmark_list[i] = iscas85_circuit_list[i];
#ifdef DEBUG_EN
				cout << selected_benchmark_list[i] << endl;
#endif
			}
			num_benchmakrs = ISCAS85_CIRCUITS;
			break;
		case CIRCUIT_GROUPS_BOOM_TILELEVEL:
			for(i=0;i<BOOM_CIRCUITS;i++) {
				selected_benchmark_list[i] = boom_circuit_list[i];
#ifdef DEBUG_EN
				cout << selected_benchmark_list[i] << endl;
#endif
			}
			num_benchmakrs = BOOM_CIRCUITS;
			break;
		default:
			cout << "Undefine index" << endl;	
	}
#ifdef ONEBENCHMARK
	num_benchmakrs = 1;
#endif
	while (infile >> b_name >> b_s >> b_c) {
		for(i=0;i<num_benchmakrs;i++) {
			if(selected_benchmark_list[i] == b_name) {
				list.Push_front(b_s, b_c, b_name, ran_round, error_rate, sqrt_win_size_min, sqrt_win_size_max);
				i = num_benchmakrs;
			}
		}
	}
#ifdef DEBUG_EN
	cout << "listnode operation test" << endl;
	cout << "Print listnode information (name,signals,cycles)" << endl;
	list.PrintList();
	cout << "Print listnode fault injection result" << endl;
	list.PrintErrors();
	cout << "Get_benchmark_listnode operation test" << endl;
	for(i=0;i<num_benchmakrs;i++) {
		pending_benchmark = list.Get_benchmark_listnode(i);
		cout << pending_benchmark.GetName() << endl;
	}
#endif
	
	cout << endl << "Evaluating..." << endl;
	//Evaluating
#ifdef MULTITHREADING
	pthread_t *thread;
	thread = (pthread_t *) malloc(sizeof(pthread_t) * (num_benchmakrs));
	benchmark_listnode *mt_pending_benchmark = new benchmark_listnode[num_benchmakrs];
	for(i=0;i<num_benchmakrs;i++) {
		mt_pending_benchmark[i] = list.Get_benchmark_listnode(i);
#ifdef DEBUG_EN
		cout << "/**********Evaluating: " << mt_pending_benchmark[i].GetName() << "**********/" << endl;
#endif
		if (pthread_create (&thread[i], NULL, MSCM,  mt_pending_benchmark+i) != 0 ) {
	      		perror("Can't create thread");
	      		free(thread);
	       		exit(-1);
	      	}
		//MSCM(pending_benchmark);
	}
	for(i=0;i<num_benchmakrs;i++)
 		pthread_join (thread[i], NULL);
#else
	for(i=0;i<num_benchmakrs;i++) {
		pending_benchmark = list.Get_benchmark_listnode(i);
		MSCM(pending_benchmark);
#ifdef DEBUG_EN
		cout << "/**********Evaluating: " << pending_benchmark.GetName() << "**********/" << endl;
#endif

	}
#endif
	

	cout << "All benchmark evaluated!!!" << endl;
	cout << "Analyze..." << endl;
	
	ResultRecordCreation(list, num_windows, num_benchmakrs);
	
	cout << endl << "Decision mode" << endl;
	if(selected_mode) { //Decision mode
		//Chooise one configuration
		pending_benchmark = list.Get_benchmark_listnode(0);
		int sqrt_win_size_index = 0;
		int win_size_number;
		int dimension_number;
#ifndef DECISION_MODE_USER_SPECIFIED
		double WeightedSumMax = 0;
		double WeightedSumMaxTemp = 0;
		int WeightedSumMaxIndex = 0;
		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {

			WeightedSumMaxTemp = 
			(pending_benchmark.result_record_list_for_each_iteration[j].Get_fcd_avg()*DECISION_MODE_W1) + 
			( (1 - pending_benchmark.result_record_list_for_each_iteration[j].Get_opa_avg()) *DECISION_MODE_W2) + 
			(pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg()*DECISION_MODE_W3);
			if((WeightedSumMaxTemp>WeightedSumMax)&&(pending_benchmark.result_record_list_for_each_iteration[j].Get_gc_avg()<DECISION_MODE_GC_BUDGET) ) {
				WeightedSumMax = WeightedSumMaxTemp;
				WeightedSumMaxIndex = j;
			}
		}
		win_size_number = sqrt( pending_benchmark.result_record_list_for_each_iteration[WeightedSumMaxIndex].GetWinSize() );
		dimension_number = pending_benchmark.result_record_list_for_each_iteration[WeightedSumMaxIndex].GetDimension();
#else
		win_size_number = DECISION_MODE_USER_SPECIFIED_SWS;
		dimension_number = DECISION_MODE_USER_SPECIFIED_DIM;
#endif
		//Generating configuration, golden signature, signature visit table
		MSCM_for_decision_mode(pending_benchmark,win_size_number,dimension_number);
	}
	list.Clear();
	delete sqrt_win_array;
	return 0;
}

void MSCM_for_decision_mode(benchmark_listnode pending_benchmark,int win_size_number,short int dimension_number) {
	string VisitPath[5] = {"Row","Col","Diagonal","Reverse Diagonal","Random"};
	short int signals = pending_benchmark.GetSignals();
	int cycles = pending_benchmark.GetCycles();
	short int num_cycle_errors = pending_benchmark.GetNumError(0);
	int i,j,k;
	int total_win;
	int cycles_expend = ceil( ((double) cycles / (double) (win_size_number*win_size_number)) ) * (win_size_number*win_size_number);
#ifndef MISR
	STAGE_ONE_COMPACTION_DATA_SIZE = 1;
	STAGE_TWO_COMPACTION_DATA_SIZE = 1;
#else
	STAGE_ONE_COMPACTION_DATA_SIZE = signals;
	STAGE_TWO_COMPACTION_DATA_SIZE = signals;
#endif
	total_win = ceil( ((double) cycles / (double) (win_size_number*win_size_number)) );
	cout << "Chosen configuration" << endl;
	cout << pending_benchmark.GetName() << endl;
	cout << cycles << endl;
	cout << signals << endl;
	cout << win_size_number << endl;
	cout << dimension_number << endl;
	string decision_mode_path = "./DecisionModeOutput/";
	string configuration_file_name = decision_mode_path + "configuration.txt";
	ofstream outfile(configuration_file_name.c_str());
	//Output configuration data
	//outfile << pending_benchmark.GetName() << endl;
	outfile << endl;
	outfile << endl;
	outfile << endl;
	outfile << cycles << endl;
	outfile << signals << endl;
	outfile << win_size_number*win_size_number << endl;
	outfile << dimension_number << endl;
	outfile << "D" << dimension_number <<" 1" << endl;
	outfile << total_win << endl;
	outfile << win_size_number << endl;
	outfile << endl;
	
	outfile << ceil( log2(total_win) ) << endl;
	outfile << ceil( log2(dimension_number) ) << endl;
	outfile << ceil( log2(cycles_expend) ) << endl;
	outfile << ceil( log2(win_size_number*win_size_number) ) << endl;
	outfile << ceil( log2(win_size_number) ) << endl;
	outfile << endl;
	unsigned lfsr_taps;
	unsigned lfsr_taps_temp;
	switch (win_size_number) {
		case 4:
			lfsr_taps = LFSR_D4;
			break;
		case 8:
			lfsr_taps = LFSR_D6;
			break;
		case 16:
			lfsr_taps = LFSR_D8;
			break;
		case 32:
			lfsr_taps = LFSR_D10;
			break;
		case 64:
			lfsr_taps = LFSR_D12;
			break;
		case 128:
			lfsr_taps = LFSR_D14;
			break;
		case 256:
			lfsr_taps = LFSR_D16;
			break;
	}
	outfile << "`LOG2_WS'b";
	for(i=win_size_number*win_size_number/2;i>0;i/=2) {
		if(i!=0)
			lfsr_taps_temp = lfsr_taps / i;  	
		outfile << lfsr_taps_temp;
		if(i!=0)		
			lfsr_taps = lfsr_taps % i;
	}
	outfile << endl;

	

	error_record **errors = pending_benchmark.errors;

	


	//Stage 1 Compaction data with fault injected
	short int **StageOneCompactionDataFI = new short int*[cycles_expend];
	short int **StageOneCompactionData = new short int*[cycles_expend];
	short int *TestDataPendingFI = new short int[signals];
	short int *TestDataPending = new short int[signals];


	for(i=0;i<cycles_expend;i++) {
		StageOneCompactionDataFI[i] = new short int[STAGE_ONE_COMPACTION_DATA_SIZE];
		StageOneCompactionData[i] = new short int[STAGE_ONE_COMPACTION_DATA_SIZE];
		for(j=0;j<STAGE_ONE_COMPACTION_DATA_SIZE;j++) {
			StageOneCompactionDataFI[i][j] = 0;
			StageOneCompactionData[i][j] = 0;
		}
	}

	//result_record temp;
	//Stage 2 Compaction data 
	short int ***StageTwoCompactionDataFI = new short int**[dimension_number];
	short int ***StageTwoCompactionData = new short int**[dimension_number];
	int ***VisitTable = new int**[dimension_number];
	short int **VisitTableCycleMapToSig = new short int*[dimension_number];
	StageTwoCompactionDataSet *GoldenSet = new StageTwoCompactionDataSet [dimension_number];
	StageTwoCompactionDataSet *FISet = new StageTwoCompactionDataSet [dimension_number];
	VisitTableSet *VisitTableSetTemp = new VisitTableSet [dimension_number];
	bool **ComparisonResults = new bool*[dimension_number];
	bool **PossibleBuggyEvent = new bool*[win_size_number];
	
	for(i=0;i<win_size_number;i++)
		PossibleBuggyEvent[i] = new bool [win_size_number];		

	for (i=0;i<dimension_number;i++) {
		StageTwoCompactionData[i] = new short int*[win_size_number];
		StageTwoCompactionDataFI[i] = new short int*[win_size_number];
		VisitTable[i] = new int*[win_size_number];
		ComparisonResults[i] = new bool [win_size_number];
		VisitTableCycleMapToSig[i] = new short int[win_size_number*win_size_number];
		for (j=0;j<win_size_number;j++) {
			StageTwoCompactionData[i][j] = new short int[STAGE_TWO_COMPACTION_DATA_SIZE];
			StageTwoCompactionDataFI[i][j] = new short int[STAGE_TWO_COMPACTION_DATA_SIZE];
			VisitTable[i][j] = new int[win_size_number];
		}
		GoldenSet[i].StageTwoCompactionData = StageTwoCompactionData[i];
		FISet[i].StageTwoCompactionData = StageTwoCompactionDataFI[i];
		VisitTableSetTemp[i].VisitTableData = VisitTable[i];
		VisitTableSetTemp[i].VisitTableCycleMapToSigData = VisitTableCycleMapToSig[i];
	}

	string circuit_name = "./signal_trace_data/signal_trace_data_" + pending_benchmark.GetName() + ".txt";
	ifstream infile(circuit_name.c_str());
	string signal_trace_data_signals;

	int ran_number = 0;
	long int win_number;
	long int cycle_number=0;
	long int error_cycles_number = 0;
	short int num_signal_errors,signal_error_number;
	int error_signals_number = 0;
	string s1_output_Golden_file_name = decision_mode_path + "s1_output_Golden.txt";
	string s1_output_FI_file_name = decision_mode_path + "s1_output_FI.txt";
	ofstream s1_output_Golden_file(s1_output_Golden_file_name.c_str());
	ofstream s1_output_FI_file(s1_output_FI_file_name.c_str());
	

	while (infile >> signal_trace_data_signals) {
		short int StageOneCompactionOUT_sim[STAGE_ONE_COMPACTION_DATA_SIZE];
		short int StageOneCompactionOUT_golden[STAGE_ONE_COMPACTION_DATA_SIZE];
		//Stage 1 compaction
		if(cycle_number!=errors[ran_number][error_cycles_number].cycle_number) {
			for(i=0;i<signals;i++) {
				//For golden
				TestDataPending[i] = signal_trace_data_signals.at(i) - '0';
			}
#ifndef MISR
			StageOneCompaction(signals,TestDataPending, StageOneCompactionOUT_golden);
#else
			for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++)
				StageOneCompactionOUT_golden[i] = TestDataPending[i];
#endif
			for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
				StageOneCompactionData[cycle_number][i] = StageOneCompactionOUT_golden[i];
				StageOneCompactionDataFI[cycle_number][i] =StageOneCompactionOUT_golden[i];
			}
		}
		else {
			num_signal_errors = errors[ran_number][error_cycles_number].list.Get_num_error();
			signal_error_number = errors[ran_number][error_cycles_number].list.Get_signal_number(error_signals_number);
			for(i=0;i<signals;i++) {
				//For golden		
				TestDataPending[i] = signal_trace_data_signals.at(i) - '0';
			
				//For sim
				TestDataPendingFI[i] = TestDataPending[i];
				if(i==signal_error_number) {
					TestDataPendingFI[i] ^= 1;
					if(error_signals_number!=(num_signal_errors-1)) {
						error_signals_number++;
						signal_error_number = errors[ran_number][error_cycles_number].list.Get_signal_number(error_signals_number);
					}
				}
			}
#ifndef MISR
			StageOneCompaction(signals, TestDataPending, StageOneCompactionOUT_golden);
			StageOneCompaction(signals, TestDataPendingFI, StageOneCompactionOUT_sim);
#else
			for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
				StageOneCompactionOUT_golden[i] = TestDataPending[i];
				StageOneCompactionOUT_sim[i] = TestDataPendingFI[i];
			}
#endif
			
			for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
				StageOneCompactionData[cycle_number][i] = StageOneCompactionOUT_golden[i];
				StageOneCompactionDataFI[cycle_number][i] =StageOneCompactionOUT_sim[i];
			}

			if(error_cycles_number!=(num_cycle_errors-1))
				error_cycles_number++;
			error_signals_number = 0;
		}
		cycle_number++;
	}
	
	//Output stage-1 compaction data
	for(i=0;i<cycles_expend;i++) {
		for(j=0;j<STAGE_ONE_COMPACTION_DATA_SIZE;j++) {
			s1_output_Golden_file << StageOneCompactionData[i][j];
			s1_output_FI_file << StageOneCompactionDataFI[i][j];	
		}
		s1_output_Golden_file << endl;
		s1_output_FI_file << endl;
	}

	total_win = ceil( ((double) cycles / (double) (win_size_number*win_size_number)) );

	VisitTableCreation(dimension_number,win_size_number,VisitTableSetTemp);	
	
	
	string sfe_output_file_name = decision_mode_path + "sfe_output.txt";
	ofstream sfe_output_file(sfe_output_file_name.c_str());

	string VisitTable_output_merge_file_name = decision_mode_path  + "Merge_VisitTable_output.txt";
	ofstream VisitTable_output_merge_file(VisitTable_output_merge_file_name.c_str());

	for(j=0;j<(win_size_number*win_size_number);j++) {
		for(i=dimension_number-1;i>=0;i--) {
			int index = VisitTableSetTemp[i].VisitTableCycleMapToSigData[j];
			for(k=win_size_number-1;k>=0;k--) {
				if(k==index)
					VisitTable_output_merge_file << "1";
				else
					VisitTable_output_merge_file << "0";
			}
		}
		VisitTable_output_merge_file << endl;
	}

	string s2_output_Golden_merge_file_name = decision_mode_path  + "Merge_s2_output_Golden.txt";
	string s2_output_FI_merge_file_name = decision_mode_path + "Merge_s2_output_FI.txt";
	string cr_output_merge_file_name = decision_mode_path + "Merge_cr_output.txt";
	string s2_in_FI_merge_file_name = decision_mode_path + "Merge_s2_in_FI.txt";

	ofstream s2_output_Golden_merge_file(s2_output_Golden_merge_file_name.c_str());
	ofstream s2_output_FI_merge_file(s2_output_FI_merge_file_name.c_str());
	ofstream cr_output_merge_file(cr_output_merge_file_name.c_str());
	ifstream s2_in_FI_merge_file(s2_in_FI_merge_file_name.c_str());
	string s2_in_FI_merge_data;
	for(win_number=0;win_number<total_win;win_number++) {
		int StageOneCompactionDataIndex = win_number * win_size_number * win_size_number;
		
		//Stage 2 compaction
		StageTwoCompaction(dimension_number,win_size_number,StageOneCompactionData+StageOneCompactionDataIndex,VisitTableSetTemp,GoldenSet);
#ifdef DECISION_MODE_CS2_DATA_FROM_FILE
		s2_in_FI_merge_file >> s2_in_FI_merge_data;
		for(i=dimension_number-1;i>=0;i--) {
			for(j=win_size_number-1;j>=0;j--) {
				for(k=0;k<STAGE_TWO_COMPACTION_DATA_SIZE;k++) {
					int index = (dimension_number*win_size_number*STAGE_TWO_COMPACTION_DATA_SIZE-1) - (i*win_size_number+j*STAGE_TWO_COMPACTION_DATA_SIZE+k);
					FISet[i].StageTwoCompactionData[j][k] = s2_in_FI_merge_data.at(index) - '0';
				}
			}
		}
#else
		StageTwoCompaction(dimension_number,win_size_number,StageOneCompactionDataFI+StageOneCompactionDataIndex,VisitTableSetTemp,FISet);
#endif
		
		
		//Comparison
		Comparison(dimension_number,win_size_number,GoldenSet,FISet,ComparisonResults);

		//Output stage-2 compaction data and compare results
		for(i=dimension_number-1;i>=0;i--) {
			for(j=win_size_number-1;j>=0;j--) {
				for(k=0;k<STAGE_TWO_COMPACTION_DATA_SIZE;k++) {
					s2_output_Golden_merge_file << GoldenSet[i].StageTwoCompactionData[j][k];
					s2_output_FI_merge_file << FISet[i].StageTwoCompactionData[j][k];
				}
				cr_output_merge_file << ComparisonResults[i][j];
			}
		}
		s2_output_Golden_merge_file << endl;
		s2_output_FI_merge_file << endl;
		cr_output_merge_file << endl;

		//Intersection
		Intersection(dimension_number,win_size_number,VisitTableSetTemp,ComparisonResults,PossibleBuggyEvent);
		
		//Output SFEs
		for(i=0;i<win_size_number;i++) {
			for(j=0;j<win_size_number;j++) {
				sfe_output_file << PossibleBuggyEvent[i][j];
				sfe_output_file << endl;
			}
		}
	}
}


#ifdef MULTITHREADING
void *MSCM(void *p_benchmark) {	
	benchmark_listnode pending_benchmark = *reinterpret_cast<benchmark_listnode *> (p_benchmark);
#else
void MSCM(benchmark_listnode pending_benchmark) {	
#endif
	string VisitPath[5] = {"Row","Col","Diagonal","Reverse Diagonal","Random"};
	short int signals = pending_benchmark.GetSignals();	
	int cycles = pending_benchmark.GetCycles();
#ifndef MISR
	STAGE_ONE_COMPACTION_DATA_SIZE = 1;
	STAGE_TWO_COMPACTION_DATA_SIZE = 1;
#else
	STAGE_ONE_COMPACTION_DATA_SIZE = signals;
	STAGE_TWO_COMPACTION_DATA_SIZE = signals;
#endif
	int sqrt_win_size_min = pending_benchmark.GetWinSizeMin();
	int sqrt_win_size_max = pending_benchmark.GetWinSizeMax();

	int cycles_expend = ceil( ((double) cycles / (double) (sqrt_win_size_max*sqrt_win_size_max)) ) * (sqrt_win_size_max*sqrt_win_size_max);
	
	int i,j,k,l;
	int ran_number=0, cycle_number=0, error_cycles_number=0, error_signals_number=0;	
	short int num_signal_errors, signal_error_number;
	int total_iteration, iteration_number;
	int win_number;
	short int dimension_number;
	string circuit_name = "./signal_trace_data/signal_trace_data_" + pending_benchmark.GetName() + ".txt";
	ifstream infile(circuit_name.c_str());
	string signal_trace_data_signals;
	error_record **errors = pending_benchmark.errors;
	int win_size_number,total_win;
	//Stage 1 Compaction data with fault injected
	short int **StageOneCompactionDataFI = new short int*[cycles_expend];
	short int **StageOneCompactionData = new short int*[cycles_expend];
	short int *TestDataPendingFI = new short int[signals];
	short int *TestDataPending = new short int[signals];

	for(i=0;i<cycles_expend;i++) {
		StageOneCompactionDataFI[i] = new short int[STAGE_ONE_COMPACTION_DATA_SIZE];
		StageOneCompactionData[i] = new short int[STAGE_ONE_COMPACTION_DATA_SIZE];
		for(j=0;j<STAGE_ONE_COMPACTION_DATA_SIZE;j++) {
			StageOneCompactionDataFI[i][j] = 0;
			StageOneCompactionData[i][j] = 0;
		}
	}

	//result_record temp;
	//Stage 2 Compaction data 
	short int ***StageTwoCompactionDataFI = new short int**[NUMBER_OF_DIMENSIONS];
	short int ***StageTwoCompactionData = new short int**[NUMBER_OF_DIMENSIONS];
	int ***VisitTable = new int**[NUMBER_OF_DIMENSIONS];
	short int **VisitTableCycleMapToSig = new short int*[NUMBER_OF_DIMENSIONS];
	StageTwoCompactionDataSet *GoldenSet = new StageTwoCompactionDataSet [NUMBER_OF_DIMENSIONS];
	StageTwoCompactionDataSet *FISet = new StageTwoCompactionDataSet [NUMBER_OF_DIMENSIONS];
	VisitTableSet *VisitTableSetTemp = new VisitTableSet [NUMBER_OF_DIMENSIONS];
	bool **ComparisonResults = new bool*[NUMBER_OF_DIMENSIONS];
	bool **PossibleBuggyEvent = new bool*[sqrt_win_size_max];
	
	for(i=0;i<sqrt_win_size_max;i++)
		PossibleBuggyEvent[i] = new bool [sqrt_win_size_max];		

	for (i=0;i<NUMBER_OF_DIMENSIONS;i++) {
		StageTwoCompactionData[i] = new short int*[sqrt_win_size_max];
		StageTwoCompactionDataFI[i] = new short int*[sqrt_win_size_max];
		VisitTable[i] = new int*[sqrt_win_size_max];
		ComparisonResults[i] = new bool [sqrt_win_size_max];
		VisitTableCycleMapToSig[i] = new short int[sqrt_win_size_max*sqrt_win_size_max];
		for (j=0;j<sqrt_win_size_max;j++) {
			StageTwoCompactionData[i][j] = new short int[STAGE_TWO_COMPACTION_DATA_SIZE];
			StageTwoCompactionDataFI[i][j] = new short int[STAGE_TWO_COMPACTION_DATA_SIZE];
			VisitTable[i][j] = new int[sqrt_win_size_max];
		}
		GoldenSet[i].StageTwoCompactionData = StageTwoCompactionData[i];
		FISet[i].StageTwoCompactionData = StageTwoCompactionDataFI[i];
		VisitTableSetTemp[i].VisitTableData = VisitTable[i];
		VisitTableSetTemp[i].VisitTableCycleMapToSigData = VisitTableCycleMapToSig[i];
	}


	for(ran_number=0;ran_number<pending_benchmark.GetRanTimes();ran_number++) {
		short int num_cycle_errors = pending_benchmark.GetNumError(ran_number);
		ifstream infile(circuit_name.c_str());
		//Stage 1 compaction
#ifdef WRITE_FAULT_INJECTED_RESULT
		string circuit_name_fault_injected = "./signal_trace_data/signal_trace_data_" + pending_benchmark.GetName() + "_FI_" + to_string(ran_number) + ".txt";
		ofstream outfile(circuit_name_fault_injected.c_str());
#endif
		cycle_number=0;
		error_cycles_number = 0;
		error_signals_number = 0;
		while (infile >> signal_trace_data_signals) {
			short int StageOneCompactionOUT_sim[STAGE_ONE_COMPACTION_DATA_SIZE];
			short int StageOneCompactionOUT_golden[STAGE_ONE_COMPACTION_DATA_SIZE];
			//Stage 1 compaction
			if(cycle_number!=errors[ran_number][error_cycles_number].cycle_number) {
				for(i=0;i<signals;i++) {
					//For golden
					TestDataPending[i] = signal_trace_data_signals.at(i) - '0';
				}
#ifndef MISR
				StageOneCompaction(signals,TestDataPending, StageOneCompactionOUT_golden);
#else
				for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++)
					StageOneCompactionOUT_golden[i] = TestDataPending[i];
#endif
#ifdef WRITE_FAULT_INJECTED_RESULT
				outfile << signal_trace_data_signals << endl;
#endif
				for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
					StageOneCompactionData[cycle_number][i] = StageOneCompactionOUT_golden[i];
					StageOneCompactionDataFI[cycle_number][i] =StageOneCompactionOUT_golden[i];
				}
			}
			else {
				num_signal_errors = errors[ran_number][error_cycles_number].list.Get_num_error();
				signal_error_number = errors[ran_number][error_cycles_number].list.Get_signal_number(error_signals_number);
				for(i=0;i<signals;i++) {
					//For golden		
					TestDataPending[i] = signal_trace_data_signals.at(i) - '0';
				
					//For sim
					TestDataPendingFI[i] = TestDataPending[i];
					if(i==signal_error_number) {
						TestDataPendingFI[i] ^= 1;
						if(error_signals_number!=(num_signal_errors-1)) {
							error_signals_number++;
							signal_error_number = errors[ran_number][error_cycles_number].list.Get_signal_number(error_signals_number);
						}
					}
#ifdef WRITE_FAULT_INJECTED_RESULT
					outfile << TestDataPendingFI[i];
#endif
				}
#ifdef WRITE_FAULT_INJECTED_RESULT
				outfile << endl;
#endif

#ifndef MISR
				StageOneCompaction(signals, TestDataPending, StageOneCompactionOUT_golden);
				StageOneCompaction(signals, TestDataPendingFI, StageOneCompactionOUT_sim);
#else
				for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
					StageOneCompactionOUT_golden[i] = TestDataPending[i];
					StageOneCompactionOUT_sim[i] = TestDataPendingFI[i];
				}
#endif
				
				for(i=0;i<STAGE_ONE_COMPACTION_DATA_SIZE;i++) {
					StageOneCompactionData[cycle_number][i] = StageOneCompactionOUT_golden[i];
					StageOneCompactionDataFI[cycle_number][i] =StageOneCompactionOUT_sim[i];
				}

				if(error_cycles_number!=(num_cycle_errors-1))
					error_cycles_number++;
				error_signals_number = 0;
			}
			cycle_number++;
		}

#ifdef DEBUG_EN
		//Fault injection & Stage 1 compaction Checking
		pending_benchmark.DisplayErrors(ran_number);
		cout << "Cycle\tGolden\tSim" << endl;
		for(i=0;i<cycles;i++)
			if(StageOneCompactionData[i][0]!=StageOneCompactionDataFI[i][0])
				cout << i << "\t" << StageOneCompactionData[i][0] << "\t" << StageOneCompactionDataFI[i][0] << endl;
#endif	
		for(win_size_number=sqrt_win_size_min ;win_size_number<=sqrt_win_size_max; win_size_number*=2) {
			total_win = ceil( ((double) cycles / (double) (win_size_number*win_size_number)) );
			total_iteration = total_win * NUMBER_OF_DIMENSIONS;

			//VisitTableCreation
			VisitTableCreation(NUMBER_OF_DIMENSIONS,win_size_number,VisitTableSetTemp);

			for(dimension_number=1;dimension_number<=NUMBER_OF_DIMENSIONS;dimension_number++) {
				long int CapturedFaultEvents=0, PossibleFaultEvents=0;
				int s1_num_aliasing=0,s2_num_aliasing=0,s2_num=0,overall_num_aliasing=0;

				for(win_number=0;win_number<total_win;win_number++) {
					int StageOneCompactionDataIndex = win_number * win_size_number * win_size_number;

					//Stage 2 compaction
					StageTwoCompaction(dimension_number,win_size_number,StageOneCompactionData+StageOneCompactionDataIndex,VisitTableSetTemp,GoldenSet);
					StageTwoCompaction(dimension_number,win_size_number,StageOneCompactionDataFI+StageOneCompactionDataIndex,VisitTableSetTemp,FISet);
					
					//Comparison
					Comparison(dimension_number,win_size_number,GoldenSet,FISet,ComparisonResults);

					//Intersection
					Intersection(dimension_number,win_size_number,VisitTableSetTemp,ComparisonResults,PossibleBuggyEvent);

					//Analyze

#ifdef ITERATION_DEBUG_EN
					int StageOneCompactionDataIndexEnd = StageOneCompactionDataIndex + (win_size_number * win_size_number) -1;
					cout << "********************This iteration information********************" << endl;
					cout << "Circuit name: " << pending_benchmark.GetName() << endl;
					cout << "Cycles: " << pending_benchmark.GetCycles() << endl;
					cout << "Signals: " << pending_benchmark.GetSignals() << endl;
					cout << "ran_number: " << ran_number << endl;
					cout << "dimension_number: " << dimension_number << endl;
					cout << "win_size_number: " << win_size_number << endl;
					cout << "win_number: " << win_number << endl;

					cout << "**********index**********" << endl;
					cout << "cn:" << pending_benchmark.GetName() << ",";
					cout << "rn:" << ran_number << ",";
					cout << "dn:" << dimension_number << ",";
					cout << "wsn:" << win_size_number << ",";
					cout << "wn:" << win_number << endl;
					cout << "**********Pending cycles**********" << endl;
					cout << "Cycle " << StageOneCompactionDataIndex << "-" << StageOneCompactionDataIndexEnd << endl;

					
					//cout << "**********Fault injection rusult after Stage 1 Compaction**********" << endl;
					//cout << "**********(Comparison of Golden and Simulation)**********" << endl;
					//cout << "Cycle\tGolden\tSim" << endl;
					//for(i=StageOneCompactionDataIndex;i<StageOneCompactionDataIndexEnd;i++)
						//if(StageOneCompactionData[i][0]!=StageOneCompactionDataFI[i][0])
							//cout << i << "\t" << StageOneCompactionData[i][0] << "\t" << StageOneCompactionDataFI[i][0] << endl;
					/*cout << "**********Stage 2 compaction visit table**********" << endl;
					for(i=1;i<=dimension_number;i++) {
						cout << VisitPath[i-1] << "path" << endl;
						for(j=0;j<win_size_number;j++) {
							cout << "S[" << j << "]: " ;
							for (k=0;k<win_size_number;k++) {
								cout << VisitTableSetTemp[i-1].VisitTableData[j][k] << " ";
							}
							cout << endl;
						}
					}*/
					
					cout << "**********Possible Buggy Event**********" << endl;
					for(i=0;i<win_size_number;i++) {
						cout << StageOneCompactionDataIndex + i*win_size_number << "-" << StageOneCompactionDataIndex + (i+1)*win_size_number -1 << ": \t";
						for(j=0;j<win_size_number;j++) {
							cout << PossibleBuggyEvent[i][j];
						}
						cout << endl;
					}
#endif
					//Analyze record
					for(i=0;i<win_size_number;i++) {
						for(j=0;j<win_size_number;j++) {
							PossibleFaultEvents += PossibleBuggyEvent[i][j];
						}
					}
					short int sig_temp[dimension_number][win_size_number];
					for(i=0;i<dimension_number;i++)
						for(j=0;j<win_size_number;j++)
							sig_temp[i][j] = 0;
					
					for(i=0;i<pending_benchmark.GetNumError(ran_number);i++) {						
						if( (errors[ran_number][i].cycle_number >= StageOneCompactionDataIndex) && (errors[ran_number][i].cycle_number <= StageOneCompactionDataIndexEnd) ) {
							k = errors[ran_number][i].cycle_number % (win_size_number*win_size_number);
							CapturedFaultEvents += PossibleBuggyEvent[k/win_size_number][k%win_size_number];

							for(j=0;j<dimension_number;j++) {
								int index = errors[ran_number][i].cycle_number % (win_size_number*win_size_number);
								k = VisitTableSetTemp[j].VisitTableCycleMapToSigData[index];
								index = errors[ran_number][i].cycle_number;
								
								for(l=0;l<STAGE_ONE_COMPACTION_DATA_SIZE;l++) {
									if(StageOneCompactionData[index][l]!=StageOneCompactionDataFI[index][l]) {
										sig_temp[j][k] = 1;
										l = STAGE_ONE_COMPACTION_DATA_SIZE;
									}
								}
								//if(StageOneCompactionData[index][0]!=StageOneCompactionDataFI[index][0])
									//sig_temp[j][k] = 1;
								//if(GoldenSet[j].StageTwoCompactionData[k][0] == FISet[j].StageTwoCompactionData[k][0])
									//s2_num_aliasing+= 1;
							}
						}
					}
					for(i=0;i<dimension_number;i++) {
						for(j=0;j<win_size_number;j++) {
							if(sig_temp[i][j]==1) {
								s2_num += 1;
								short int s1_non_aliasing_flag = 0;
								for(k=0;k<STAGE_TWO_COMPACTION_DATA_SIZE;k++) {
									if(GoldenSet[i].StageTwoCompactionData[j][k] != FISet[i].StageTwoCompactionData[j][k]) {
										s1_non_aliasing_flag = 1;
									}
								}
								if(s1_non_aliasing_flag==0) {
									cout << "s2 Aliasing occured!!!" << endl;
									cout << "D" << i << ", " << "Sig_num" << j << endl;
									cout << "wn:" << win_number << endl;
									s2_num_aliasing+= 1;
								}
							}
						}
					}

					//break; //First window
				}
				//break; // One dimension

				//Analyze record
				double FaultCycleDensity;
				double s1_pa, s2_pa, o_pa;
				double cr;
				float gc=0;
				if(PossibleFaultEvents!=0)
					FaultCycleDensity = (double) CapturedFaultEvents / PossibleFaultEvents;
				else {
					FaultCycleDensity = (double) 1;
					cout << "FaultCycleDensity Bug occur!!!" << endl;
					cout << "PossibleFaultEvents = 0" << endl;
				}
				if(PossibleFaultEvents<CapturedFaultEvents) {
					cout << "FaultCycleDensity Bug occur!!!" << endl;
					cout << "PossibleFaultEvents < CapturedFaultEvents" << endl;
				}
				for(i=0;i<pending_benchmark.GetNumError(ran_number);i++) {
					k = errors[ran_number][i].cycle_number;
					short int s1_non_aliasing_flag = 0;
					for(j=0;j<STAGE_ONE_COMPACTION_DATA_SIZE;j++) {
						if(StageOneCompactionData[k][j]!=StageOneCompactionDataFI[k][j]) {
							s1_non_aliasing_flag = 1;
						}
					}
					if(s1_non_aliasing_flag==0) {
						s1_num_aliasing+= 1;
					}
				}
				s1_pa = (double) s1_num_aliasing / pending_benchmark.GetNumError(ran_number);
				s2_pa = (double) s2_num_aliasing / s2_num;
				o_pa = (double) CapturedFaultEvents / pending_benchmark.GetNumError(ran_number);
				o_pa = (double) (1 - o_pa);
				cr = (double) ( (double) (win_size_number * dimension_number * total_win) / (double) (signals*cycles) );
				cr = (double) (1 - cr);
				if(win_size_number==16)
					gc = 737.1837533 + (signals * 26.99796529) + (cycles * 0.006250714) + (dimension_number * 1509.8);
				else if(win_size_number==32)
					gc = 756.8125401 + (signals * 26.9633542) + (cycles * 0.002142857) + (dimension_number * 2701.85);
				else if(win_size_number==64)
					gc = 334.5426256 + (signals * 27.13588801) + (cycles * 0.0137975) + (dimension_number * 5037.158333);
				gc /= STANDARD_CELL_BASE_AREA;
				i = ((log2(win_size_number) - log2(sqrt_win_size_min))*NUMBER_OF_DIMENSIONS) + (dimension_number-1);
				pending_benchmark.result_record_list_for_each_iteration[i].Push_front(ran_number,FaultCycleDensity,cr,s1_pa,s2_pa,o_pa,gc);
#ifdef ITERATION_DEBUG_EN
				cout << "********************Analyze record********************" << endl;
				cout << "********************This iteration information********************" << endl;
				cout << "Circuit name: " << pending_benchmark.GetName() << endl;
				cout << "Cycles: " << pending_benchmark.GetCycles() << endl;
				cout << "Signals: " << pending_benchmark.GetSignals() << endl;
				cout << "ran_number: " << ran_number << endl;
				cout << "dimension_number: " << dimension_number << endl;
				cout << "win_size_number: " << win_size_number << endl;
				cout << "win_number: " << win_number << endl;
				cout << "**********index**********" << endl;
				cout << "cn:" << pending_benchmark.GetName() << ",";
				cout << "rn:" << ran_number << ",";
				cout << "dn:" << dimension_number << ",";
				cout << "wsn:" << win_size_number << ",";
				cout << "wn:" << win_number << endl;
				cout << "**********Analyze result**********" << endl;
				cout << "Possible Fault Events: " << PossibleFaultEvents << endl;
				cout << "Captured Fault Events: " << CapturedFaultEvents << endl;
				cout << "Fault Cycle Density: " << FaultCycleDensity << endl;
				cout << "S1_aliasing_count: " << s1_num_aliasing << endl;
				cout << "S1_pa: " << s1_pa << endl;
				cout << "S2_aliasing_count: " << s2_num_aliasing << endl;
				cout << "s2_num: " << s2_num << endl;
				cout << "S2_pa: " << s2_pa << endl;
				cout << "Overall_pa: " << o_pa << endl;
				cout << "Compression ratio: " << cr << endl;
				cout << "Gate conut: " << gc << endl;
				cout << "**********Fault list**********" << endl;
				pending_benchmark.DisplayErrors(ran_number);
#endif
			}
			//break; //win_size_small
		}

#ifdef DEBUG_EN
		cout << "ran " << ran_number << " finished" << endl << endl;
#endif
		//break; //
	}
	delete StageOneCompactionData;
	delete StageOneCompactionDataFI;
	delete TestDataPending;
	delete TestDataPendingFI;
	delete StageTwoCompactionDataFI;
	delete StageTwoCompactionData;
	delete VisitTable;

	delete GoldenSet;
	delete FISet;
	delete VisitTableSetTemp;
	delete PossibleBuggyEvent;
	//return temp;
}

void ResultRecordCreation(benchmark_linkedlist list,int num_windows, int num_benchmakrs) {
	double avg_fcd[num_windows*NUMBER_OF_DIMENSIONS];
	double avg_s1pa[num_windows*NUMBER_OF_DIMENSIONS];
	double avg_s2pa[num_windows*NUMBER_OF_DIMENSIONS];
	double avg_opa[num_windows*NUMBER_OF_DIMENSIONS];
	double avg_cr[num_windows*NUMBER_OF_DIMENSIONS];
	int win_size_record[num_windows*NUMBER_OF_DIMENSIONS];
	int dimension_record[num_windows*NUMBER_OF_DIMENSIONS];
	benchmark_listnode pending_benchmark;
	int i,j;

	//Initialize
	pending_benchmark = list.Get_benchmark_listnode(0);
	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		cout << "win_size: " << pending_benchmark.result_record_list_for_each_iteration[j].GetWinSize() << endl;
		win_size_record[j] = pending_benchmark.result_record_list_for_each_iteration[j].GetWinSize();
		cout << "dimension: " << pending_benchmark.result_record_list_for_each_iteration[j].GetDimension() << endl;
		dimension_record[j] = pending_benchmark.result_record_list_for_each_iteration[j].GetDimension();
		avg_fcd[j] = 0;
		avg_s1pa[j] = 0;
		avg_s2pa[j] = 0;
		avg_opa[j] = 0;
		avg_cr[j] = 0;
	}

	for(i=0;i<num_benchmakrs;i++) {
		pending_benchmark = list.Get_benchmark_listnode(i);
		cout << "----------------------" << pending_benchmark.GetName() << "---------------------" << endl;
#ifdef VOTING
		string benchmark_result_record_name = "Result_records_voting_" + pending_benchmark.GetName() +".txt";
#else
		string benchmark_result_record_name = "Result_records_" + pending_benchmark.GetName() +".txt";
#endif
		string benchmark_result_reocrd_file = "./Result_records/" + benchmark_result_record_name;
		ofstream benchmark_outfile(benchmark_result_reocrd_file.c_str());
		
		cout << "Avg. Fault Cycle Density" << endl;
		benchmark_outfile << "Avg. Fault Cycle Density" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_fcd_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_fcd_avg() << "\t";
		}
		cout << endl;
		benchmark_outfile << endl;
		cout << "Avg. S1_pa" << endl;
		benchmark_outfile << "Avg. S1_pa" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_s1pa_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_s1pa_avg() << "\t";
		}

		cout << endl;
		benchmark_outfile << endl;
		cout << "Avg. S2_pa" << endl;
		benchmark_outfile << "Avg. S2_pa" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_s2pa_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_s2pa_avg() << "\t";
		}

		cout << endl;
		benchmark_outfile << endl;
		cout << "Avg. Overall_pa:" << endl;
		benchmark_outfile << "Avg. Overall_pa:" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_opa_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_opa_avg() << "\t";
		}

		cout << endl;
		benchmark_outfile << endl;
		cout << "Avg. Compression ratio:" << endl;
		benchmark_outfile << "Avg. Compression ratio:" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg() << "\t";
		}

		cout << endl;
		benchmark_outfile << endl;
		cout << "Avg. Gate count:" << endl;
		benchmark_outfile << "Avg. Gate count:" << endl;
		cout << "\t\t" ;
		benchmark_outfile << "\t\t" ;
		for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
			cout <<  j+1 << "\t\t";
			benchmark_outfile <<  j+1 << "\t\t";
		}

		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			if(j%NUMBER_OF_DIMENSIONS==0) {
				cout << endl;
				benchmark_outfile << endl;
				cout << win_size_record[j] << "\t";
				benchmark_outfile << win_size_record[j] << "\t";
			}
			cout << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_gc_avg() << "\t";
			benchmark_outfile << setw(9) << pending_benchmark.result_record_list_for_each_iteration[j].Get_gc_avg() << "\t";
		}



		for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
			cout << "---Analyze results---" << endl;
			
			cout << "Avg. Fault Cycle Density: " << pending_benchmark.result_record_list_for_each_iteration[j].Get_fcd_avg() << endl;
			avg_fcd[j] += pending_benchmark.result_record_list_for_each_iteration[j].Get_fcd_avg();
			cout << "Avg. S1_pa: " << pending_benchmark.result_record_list_for_each_iteration[j].Get_s1pa_avg() << endl;
			avg_s1pa[j] += pending_benchmark.result_record_list_for_each_iteration[j].Get_s1pa_avg();
			cout << "Avg. S2_pa: " << pending_benchmark.result_record_list_for_each_iteration[j].Get_s2pa_avg() << endl;
			avg_s2pa[j] += pending_benchmark.result_record_list_for_each_iteration[j].Get_s2pa_avg();
			cout << "Avg. Overall_pa: " << pending_benchmark.result_record_list_for_each_iteration[j].Get_opa_avg() << endl;
			avg_opa[j] += pending_benchmark.result_record_list_for_each_iteration[j].Get_opa_avg();
			cout << "Avg. Compression ratio: " << pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg() << endl;
			avg_cr[j] += pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg();
			if(pending_benchmark.result_record_list_for_each_iteration[j].Get_cr_avg()>1)
				cout << "Fault flag!!!!!!" << endl;
			cout << "avg_cr[" << j << "]: " << avg_cr[j] << endl;

		}
		
	}
#ifdef VOTING
	string result_record_name = "Result_records_voting.txt";
#else
	string result_record_name = "Result_records.txt";
#endif
	string result_reocrd_file = "./Result_records/" + result_record_name;
	ofstream outfile(result_reocrd_file.c_str());
	cout << "---All benchmakr analyze results---" << endl;
	cout << "num_benchmakrs: " << num_benchmakrs << endl;
	
	cout << "Avg. Fault Cycle Density" << endl;
	outfile << "Avg. Fault Cycle Density" << endl;
	cout << "\t\t" ;
	outfile << "\t\t" ;
	for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
		cout <<  j+1 << "\t\t";
		outfile <<  j+1 << "\t\t";
	}

	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		if(j%NUMBER_OF_DIMENSIONS==0) {
			cout << endl;
			outfile << endl;
			cout << win_size_record[j] << "\t";
			outfile << win_size_record[j] << "\t";
		}
		avg_fcd[j] /= num_benchmakrs;
		cout << setw(9) << avg_fcd[j] << "\t";
		outfile << setw(9) << avg_fcd[j] << "\t";
	}
	cout << endl;
	outfile << endl;
	cout << "Avg. S1_pa" << endl;
	outfile << "Avg. S1_pa" << endl;
	cout << "\t\t" ;
	outfile << "\t\t" ;
	for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
		cout <<  j+1 << "\t\t";
		outfile <<  j+1 << "\t\t";
	}

	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		if(j%NUMBER_OF_DIMENSIONS==0) {
			cout << endl;
			outfile << endl;
			cout << win_size_record[j] << "\t";
			outfile << win_size_record[j] << "\t";
		}
		avg_s1pa[j] /= num_benchmakrs;
		cout << setw(9) << avg_s1pa[j] << "\t";
		outfile << setw(9) << avg_s1pa[j] << "\t";
	}

	cout << endl;
	outfile << endl;
	cout << "Avg. S2_pa" << endl;
	outfile << "Avg. S2_pa" << endl;
	cout << "\t\t" ;
	outfile << "\t\t" ;
	for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
		cout <<  j+1 << "\t\t";
		outfile <<  j+1 << "\t\t";
	}

	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		if(j%NUMBER_OF_DIMENSIONS==0) {
			cout << endl;
			outfile << endl;
			cout << win_size_record[j] << "\t";
			outfile << win_size_record[j] << "\t";
		}
		avg_s2pa[j] /= num_benchmakrs;
		cout << setw(9) << avg_s2pa[j] << "\t";
		outfile << setw(9) << avg_s2pa[j] << "\t";
	}

	cout << endl;
	outfile << endl;
	cout << "Avg. Overall_pa:" << endl;
	outfile << "Avg. Overall_pa:" << endl;
	cout << "\t\t" ;
	outfile << "\t\t" ;
	for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
		cout <<  j+1 << "\t\t";
		outfile <<  j+1 << "\t\t";
	}

	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		if(j%NUMBER_OF_DIMENSIONS==0) {
			cout << endl;
			outfile << endl;
			cout << win_size_record[j] << "\t";
			outfile << win_size_record[j] << "\t";
		}
		avg_opa[j] /= num_benchmakrs;
		cout << setw(9) << avg_opa[j] << "\t";
		outfile << setw(9) << avg_opa[j] << "\t";
	}

	cout << endl;
	outfile << endl;
	cout << "Avg. Compression ratio:" << endl;
	outfile << "Avg. Compression ratio:" << endl;
	cout << "\t\t" ;
	outfile << "\t\t" ;
	for(j=0;j<(NUMBER_OF_DIMENSIONS);j++) {
		cout <<  j+1 << "\t\t";
		outfile <<  j+1 << "\t\t";
	}

	for(j=0;j<(num_windows*NUMBER_OF_DIMENSIONS);j++) {
		if(j%NUMBER_OF_DIMENSIONS==0) {
			cout << endl;
			outfile << endl;
			cout << win_size_record[j] << "\t";
			outfile << win_size_record[j] << "\t";
		}
		avg_cr[j] /= num_benchmakrs;
		cout << setw(9) << avg_cr[j] << "\t";
		outfile << setw(9) << avg_cr[j] << "\t";
	}
}

void StageOneCompaction(short int Signals,short int* signal_trace_data_signals,short int* output) {
	int i;	
	*output = 0;
	for(i=0;i<Signals;i++)
		*output ^= signal_trace_data_signals[i];
}

void Intersection(short int dimension_number, int sqrt_win_size, VisitTableSet *VisitTable, bool **ComparisonResults, bool **PossibleBuggyEvent) {
	int i,j,k,index,count;
	bool temp=1;
	for(i=0;i<sqrt_win_size;i++) {
		for(j=0;j<sqrt_win_size;j++) {
			temp = 1;
			count = 0;
			for(k=0;k<dimension_number;k++) {
				index = VisitTable[k].VisitTableCycleMapToSigData[(i*sqrt_win_size)+j];
				//cout << "index: " << index << endl;
				count += (ComparisonResults[k][index]==1);
				temp &= ComparisonResults[k][index];
				//cout << "ComparisonResults: " << ComparisonResults[k][index] << endl;
			}
#ifdef VOTING
			int threshole;
			/*if((dimension_number%2)==0)
				threshole = ((dimension_number/2));
			else
				threshole = ((dimension_number/2) +1);*/
			threshole = dimension_number - 4;
			if(dimension_number>4)
				temp = ( count >= threshole ? 1 : 0);
#endif
			PossibleBuggyEvent[i][j] = temp;
		}
	}

	//cout << "Intersection" << endl;
}

void Comparison(short int dimension_number,int sqrt_win_size,StageTwoCompactionDataSet *GoldenSet,StageTwoCompactionDataSet *FISet,bool **ComparisonResults) {
	//cout << "Comparison" << endl;
	int i,j,k;
	for(i=0;i<dimension_number;i++) {
		for(j=0;j<sqrt_win_size;j++) {
			//ComparisonResults[i][j] = GoldenSet[i].StageTwoCompactionData[j][0] != FISet[i].StageTwoCompactionData[j][0];
			ComparisonResults[i][j] = false;
			for(k=0;k<STAGE_TWO_COMPACTION_DATA_SIZE;k++) {
				if(GoldenSet[i].StageTwoCompactionData[j][k] != FISet[i].StageTwoCompactionData[j][k]) {
					ComparisonResults[i][j] = true;
					k = STAGE_TWO_COMPACTION_DATA_SIZE;
				}
			}
		}
	}
}

void StageTwoCompaction(short int dimension_number,int sqrt_win_size,short int **StageOneCompactionDataIN, VisitTableSet *VisitTable, StageTwoCompactionDataSet *StageTwoCompactionDataOUT) {
	int i,j;	
	//cout << "dimension_number: " << dimension_number << endl;
	//cout << "Cycle number contents" << endl;
	//for(i=0;i<sqrt_win_size;i++) {
			//for(j=0;j<sqrt_win_size;j++)
				//cout << (i*sqrt_win_size) +j << "\t";
			//cout << endl;
	//}
#ifndef MISR
	switch (dimension_number) {
		case 5:	//5D Random Visit Compaction
			StageTwoXoredCompaction(sqrt_win_size,VisitTable[4].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[4].StageTwoCompactionData);
		case 4:	//4D Reverse Diagonal Visit Compaction
			StageTwoXoredCompaction(sqrt_win_size,VisitTable[3].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[3].StageTwoCompactionData);
		case 3:	//3D Diagonal Visit Compaction
			StageTwoXoredCompaction(sqrt_win_size,VisitTable[2].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[2].StageTwoCompactionData);
		case 2:	//2D Col Visit Compaction
			StageTwoXoredCompaction(sqrt_win_size,VisitTable[1].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[1].StageTwoCompactionData);
		case 1:	//1D Row Visit Compaction
			StageTwoXoredCompaction(sqrt_win_size,VisitTable[0].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[0].StageTwoCompactionData);
			break;
		default:
			cout << "Undefined dimension!!!" << endl;
	}
#else
	switch (dimension_number) {
		case 5:	//5D Random Visit Compaction
			StageTwoMISRCompaction(sqrt_win_size,VisitTable[4].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[4].StageTwoCompactionData);
		case 4:	//4D Reverse Diagonal Visit Compaction
			StageTwoMISRCompaction(sqrt_win_size,VisitTable[3].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[3].StageTwoCompactionData);
		case 3:	//3D Diagonal Visit Compaction
			StageTwoMISRCompaction(sqrt_win_size,VisitTable[2].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[2].StageTwoCompactionData);
		case 2:	//2D Col Visit Compaction
			StageTwoMISRCompaction(sqrt_win_size,VisitTable[1].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[1].StageTwoCompactionData);
		case 1:	//1D Row Visit Compaction
			StageTwoMISRCompaction(sqrt_win_size,VisitTable[0].VisitTableData, StageOneCompactionDataIN,  StageTwoCompactionDataOUT[0].StageTwoCompactionData);
			break;
		default:
			cout << "Undefined dimension!!!" << endl;
	}
#endif
}

void VisitTableCreation(short int dimension_number, int sqrt_win_size, VisitTableSet *VisitTable) {
	switch (dimension_number) {
		case 5:	//5D Random VisitTableCreation
			//cout << "5D Random VisitTableCreation..." << endl;
			RanVisitTableCreation(sqrt_win_size, &VisitTable[4]);
		case 4:	//4D Reverse Diagonal VisitTableCreation
			//cout << "4D Reverse Diagonal VisitTableCreation.." << endl;
			RevDiaVisitTableCreation(sqrt_win_size, &VisitTable[3]);
		case 3:	//3D Diagonal VisitTableCreation
			//cout << "3D Diagonal VisitTableCreation..." << endl;
			DiaVisitTableCreation(sqrt_win_size, &VisitTable[2]);
		case 2:	//2D Col VisitTableCreation
			//cout << "2D Col VisitTableCreation..." << endl;
			ColVisitTableCreation(sqrt_win_size, &VisitTable[1]);
		case 1:	//1D Row VisitTableCreation
			//cout << "1D Row VisitTableCreation..." << endl;
			RowVisitTableCreation(sqrt_win_size, &VisitTable[0]);
			break;
		default:
			cout << "Undefined dimension!!!" << endl;
	}
}

void RowVisitTableCreation(int sqrt_win_size,VisitTableSet *VisitTable) {
	int i,j,sig_index,cycle_index;
	for(i=0;i<sqrt_win_size;i++) {
		for (j=0;j<sqrt_win_size;j++) {
			cycle_index = (i*sqrt_win_size + j);
			sig_index = i;
			VisitTable->VisitTableData[i][j] = cycle_index;
			VisitTable->VisitTableCycleMapToSigData[cycle_index] = sig_index;
		}
	}
}
void ColVisitTableCreation(int sqrt_win_size,VisitTableSet *VisitTable) {
	int i,j,sig_index,cycle_index;
	for(i=0;i<sqrt_win_size;i++) {
		for (j=0;j<sqrt_win_size;j++) {
			cycle_index = (j*sqrt_win_size + i);
			sig_index = i;
			VisitTable->VisitTableData[i][j] = cycle_index;
			VisitTable->VisitTableCycleMapToSigData[cycle_index] = sig_index;
		}
	}
}
void DiaVisitTableCreation(int sqrt_win_size,VisitTableSet *VisitTable) {
	int i,j,Diagonal_line_index,cycle_index,sig_index;
	for(i=0;i<sqrt_win_size;i++) {
		Diagonal_line_index = (i*sqrt_win_size) + i;

		for (j=0;j<sqrt_win_size;j++) {
			cycle_index = (i * sqrt_win_size) + j;
			if(cycle_index==Diagonal_line_index) { //Diagonal 
				VisitTable->VisitTableData[0][j] = cycle_index;	
				sig_index = 0;	
			}
			else if(cycle_index>Diagonal_line_index) { //Up 
				VisitTable->VisitTableData[j-i][i] = cycle_index;
				sig_index = j-i;
			}
			else { //Down 
				VisitTable->VisitTableData[sqrt_win_size-i+j][i] = cycle_index;
				sig_index = sqrt_win_size-i+j;
			}
			VisitTable->VisitTableCycleMapToSigData[cycle_index] = sig_index;
		}
	}
}
void RevDiaVisitTableCreation(int sqrt_win_size,VisitTableSet *VisitTable) {
	int i,j,Reverse_diagonal_line_index,cycle_index,sig_index;
	for(i=0;i<sqrt_win_size;i++) {
		Reverse_diagonal_line_index = (sqrt_win_size-1) * (i+1);
		for (j=0;j<sqrt_win_size;j++) {
			cycle_index = (i * sqrt_win_size) + j;
			if(cycle_index==Reverse_diagonal_line_index) { //Diagonal 
				VisitTable->VisitTableData[(sqrt_win_size-1)][i] = cycle_index;	
				sig_index = (sqrt_win_size-1);		
			}
			else if(cycle_index<Reverse_diagonal_line_index) { //Up 
				VisitTable->VisitTableData[j+i][i] = cycle_index;
				sig_index = j+i;
			}
			else { //Down 
				VisitTable->VisitTableData[(i+j)-sqrt_win_size][i] = cycle_index;
				sig_index = (i+j)-sqrt_win_size;
			}
			VisitTable->VisitTableCycleMapToSigData[cycle_index] = sig_index;
		}
	}
}
void RanVisitTableCreation(int sqrt_win_size,VisitTableSet *VisitTable) {
	int i,j,index,cycle_index,sig_index;
	int sig_count[sqrt_win_size];
	unsigned start_state = 0x1u;
	unsigned lfsr_taps;
	unsigned lfsr = start_state;
	//unsigned lsb = 0;

	switch (sqrt_win_size) {
		case 4:
			lfsr_taps = LFSR_D4;
			break;
		case 8:
			lfsr_taps = LFSR_D6;
			break;
		case 16:
			lfsr_taps = LFSR_D8;
			break;
		case 32:
			lfsr_taps = LFSR_D10;
			break;
		case 64:
			lfsr_taps = LFSR_D12;
			break;
		case 128:
			lfsr_taps = LFSR_D14;
			break;
		case 256:
			lfsr_taps = LFSR_D16;
			break;
	}

	for(i=0;i<sqrt_win_size;i++)
		sig_count[i] = 0;

	for(i=0;i<sqrt_win_size;i++) {
		for (j=0;j<sqrt_win_size;j++) {
			cycle_index = (i * sqrt_win_size) + j;
			if(cycle_index==0) {
				sig_index = 0;
				
			}	
			else {
				unsigned lsb = lfsr & 1;
				sig_index = lfsr % sqrt_win_size;
				lfsr >>= 1;
				if(lsb == 1)
					lfsr ^= lfsr_taps;
			}
			VisitTable->VisitTableData[sig_index][sig_count[sig_index]] = cycle_index;
			VisitTable->VisitTableCycleMapToSigData[cycle_index] = sig_index;
			sig_count[sig_index]++;
		}
	}
}

void StageTwoXoredCompaction(int sqrt_win_size,int **VisitTable,short int** StageOneCompactionDataIN,short int** StageTwoCompactionDataOUT) {
	int i, j, k, index;	
	short int output[STAGE_ONE_COMPACTION_DATA_SIZE];
	
	for(i=0;i<sqrt_win_size;i++) {
		for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++)
			output[k] = 0;
		for(j=0;j<sqrt_win_size;j++) {
			index = VisitTable[i][j];
			for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++)
				output[k] ^= StageOneCompactionDataIN[index][k];
		}
		for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++)
			StageTwoCompactionDataOUT[i][k] = output[k];
	}
}

void StageTwoMISRCompaction(int sqrt_win_size,int **VisitTable,short int** StageOneCompactionDataIN,short int** StageTwoCompactionDataOUT) {
	int i, j, k, index;	
	short int output[STAGE_ONE_COMPACTION_DATA_SIZE];
	short int output_temp[STAGE_ONE_COMPACTION_DATA_SIZE];
	short int MISR_taps[STAGE_ONE_COMPACTION_DATA_SIZE];
	for(i=STAGE_ONE_COMPACTION_DATA_SIZE-1;i>=0;i--) {
		if(i>=STAGE_ONE_COMPACTION_DATA_SIZE-2) {
			MISR_taps[i] = 1;		
		}	
		else {
			MISR_taps[i] = 0;
		}
	}
	for(i=0;i<sqrt_win_size;i++) {
		for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++)
			output[k] = 0;
		for(j=0;j<sqrt_win_size;j++) {
			index = VisitTable[i][j];
			for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++) {
				output_temp[k] = output[k];			
			}
			for(k=0;k<STAGE_ONE_COMPACTION_DATA_SIZE;k++) {
				if(k==STAGE_ONE_COMPACTION_DATA_SIZE-1) {
					output[k] = StageOneCompactionDataIN[index][k];
				}
				else {
					output[k] = StageOneCompactionDataIN[index][k]^output_temp[k+1];
				}
				if(MISR_taps[k]==1) 
					output[k] ^= output_temp[0];
					
			}
		}
		for(k=0;k<STAGE_TWO_COMPACTION_DATA_SIZE;k++)
			StageTwoCompactionDataOUT[i][k] = output[k];
	}
}


