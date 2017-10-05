#define ISCAS85_CIRCUITS 10
#define BOOM_CIRCUITS 1
#define CIRCUIT_GROUPS 6
#define BOOM_CIRCUITS_INDEX ISCAS85_CIRCUITS
#define CIRCUIT_GROUPS_INDEX BOOM_CIRCUITS_INDEX+BOOM_CIRCUITS
#define CIRCUIT_GROUPS_ISCAS85 ISCAS85_CIRCUITS + BOOM_CIRCUITS + 1
#define CIRCUIT_GROUPS_BOOM_TILELEVEL CIRCUIT_GROUPS_ISCAS85 + 1
#define NUMBER_OF_DIMENSIONS 5
#define LFSR_D4 0x9
#define LFSR_D6 0x2D
#define LFSR_D8 0x8E
#define LFSR_D10 0x204
#define LFSR_D12 0x829
#define LFSR_D14 0x2015
#define LFSR_D16 0x8016
#define MAXIMUM_CYCLES 100000
#define MAXIMUM_SIGNALS 1000
#define WINDOW_SPACING 5
//#define STAGE_ONE_COMPACTION_DATA_SIZE 1
//#define STAGE_TWO_COMPACTION_DATA_SIZE 1
//#define MISR 1
#define FIXED_INPUTS 1
#define NUM_THREAD 10
//#define MULTITHREADING 1
#define ONEBENCHMARK 1
//#define DEBUG_EN 1
#define ITERATION_DEBUG_EN 1
//#define FAULT_INJECTION_FROM_FILE 1
#define WRITE_FAULT_INJECTED_RESULT 1
#define VOTING 1
#define INTERSECTING_DEGREE 1
//#define DECISION_MODE_CS2_DATA_FROM_FILE 1
#define DECISION_MODE_USER_SPECIFIED 1
#define DECISION_MODE_USER_SPECIFIED_SWS 64
#define DECISION_MODE_USER_SPECIFIED_DIM 5
//#define DECISION_MODE_USER_SPECIFIED
//#define DECISION_MODE_GC_BUDGET 140000
//#define DECISION_MODE_W1 1
//#define DECISION_MODE_W2 5
//#define DECISION_MODE_W3 1
#define STANDARD_CELL_BASE_AREA 2.8224


using namespace std;


//Error_record
class error_record_for_signals_linkedlist;

class error_record_for_signals {
	private:
		int signal_number;
		error_record_for_signals *next;
	public:
		error_record_for_signals():signal_number(0) {};
		error_record_for_signals(int n):signal_number(n) {};
	friend class error_record_for_signals_linkedlist;	
};

class error_record_for_signals_linkedlist {
	private:
		error_record_for_signals *first;
		int num_error;
	public:
		error_record_for_signals_linkedlist():first(0) {};
    		void Push_front(int n);
		int Get_signal_number(int index);
		void Set_num_error(int n);
		int Get_num_error();
    		void Clear();
};

void error_record_for_signals_linkedlist::Push_front(int n){
    error_record_for_signals *newNode = new error_record_for_signals(n);   
    newNode->next = first;                 
    first = newNode;                       
}

int error_record_for_signals_linkedlist::Get_signal_number(int index){
    int i;
    error_record_for_signals IndexNode = *first;
    for(i=0;i<index;i++)
	IndexNode = *(IndexNode.next);
    return IndexNode.signal_number;                      
}

void error_record_for_signals_linkedlist::Set_num_error(int n){
    num_error = n;                    
}

int error_record_for_signals_linkedlist::Get_num_error(){
    return num_error;                    
}

void error_record_for_signals_linkedlist::Clear(){

    while (first != 0) {            // Traversal
        error_record_for_signals *current = first;
        first = first->next;
        delete current;
        current = 0;
    }
}

struct error_record {
	int cycle_number;
	error_record_for_signals_linkedlist list;
};

//Result record
class result_record_linkedlist;

class result_record {
	private:
		
		short int ran_times;
		double fault_cycle_density;
		double compression_ratio; 
		double s1_probability_of_aliasing;
		double s2_probability_of_aliasing;
		double overall_probability_of_aliasing;
		float gate_count;

		result_record *next;
	public:
		
		result_record():ran_times(0),fault_cycle_density(0),compression_ratio(0),s1_probability_of_aliasing(0),s2_probability_of_aliasing(0),overall_probability_of_aliasing(0),gate_count(0) {};
		result_record(short int rt,double fcd,double cr,double s1_pa,double s2_pa,double o_pa,float gc):ran_times(rt),fault_cycle_density(fcd),compression_ratio(cr),s1_probability_of_aliasing(s1_pa),s2_probability_of_aliasing(s2_pa),overall_probability_of_aliasing(o_pa),gate_count(gc) {};
	friend class result_record_linkedlist;	
};

class result_record_linkedlist {
	private:
		result_record *first;
		int win_size;
		short int dimension;
	public:
		result_record_linkedlist():first(NULL) {};
    		void Push_front(short int rt,double fcd,double cr,double s1_pa,double s2_pa,double o_pa,float gc);
		void SetWinSize(int sws);
		void SetDimension(short int sd);
		int GetWinSize();
		short int GetDimension();
		double Get_fcd_avg();
		double Get_cr_avg();
		double Get_s1pa_avg();
		double Get_s2pa_avg();
		double Get_opa_avg();
		double Get_gc_avg();
};

int result_record_linkedlist::GetWinSize(){
	return win_size;
}

short int result_record_linkedlist::GetDimension(){
	return dimension;
}

void result_record_linkedlist::SetWinSize(int sws){
	win_size = sws;
}

void result_record_linkedlist::SetDimension(short int sd){
	dimension = sd;
}

double result_record_linkedlist::Get_fcd_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.fault_cycle_density;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.fault_cycle_density;
		count += 1;
		
	}
	avg /= count;
	return avg;   
}

double result_record_linkedlist::Get_cr_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.compression_ratio;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.compression_ratio;
		count += 1;
		
	}
	avg /= count;
	return avg;
}

double result_record_linkedlist::Get_s1pa_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.s1_probability_of_aliasing;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.s1_probability_of_aliasing;
		count += 1;
		
	}
	avg /= count;
	return avg;   
}

double result_record_linkedlist::Get_s2pa_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.s2_probability_of_aliasing;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.s2_probability_of_aliasing;
		count += 1;
		
	}
	avg /= count;
	return avg;   
}

double result_record_linkedlist::Get_opa_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.overall_probability_of_aliasing;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.overall_probability_of_aliasing;
		count += 1;
		
	}
	avg /= count;
	return avg;   
}

double result_record_linkedlist::Get_gc_avg() {
	int i,count=0;
	double avg=0;
	result_record Node = *first;
	avg += Node.gate_count;
	count += 1;
	while((Node.next)!=NULL) {
		Node = *(Node.next);
		avg += Node.gate_count;
		count += 1;
	}
	avg /= count;
	return avg;   
}

void result_record_linkedlist::Push_front(short int rt,double fcd,double cr,double s1_pa,double s2_pa,double o_pa,float gc){
    result_record *newNode = new result_record(rt,fcd,cr,s1_pa,s2_pa,o_pa,gc);   
    newNode->next = first;                 
    first = newNode;                       
}

//Benchmark information struct
class benchmark_linkedlist;

class benchmark_listnode {
	private:
		int signals;
		int cycles;
		string name;
		benchmark_listnode *next;
		int *num_error;
		int ran_times;
		int sqrt_win_size_min;
		int sqrt_win_size_max;
	public:
		error_record **errors;
		result_record_linkedlist *result_record_list_for_each_iteration;
		benchmark_listnode():signals(0),cycles(0),name(""),next(0) {};
		benchmark_listnode(int s,int c,string n,int ran_round ,double *error_rate,int min,int max);
		int GetSignals();
		int GetCycles();
		int GetNumError(int);
		int GetRanTimes();
		int GetWinSizeMin();
		int GetWinSizeMax();
		//void Ceration_result_record(int num_win_size);
		//void Push_result_record(result_record_array);
		string GetName();
		void DisplayErrors();
		void DisplayErrors(int ran_time_number);

	friend class benchmark_linkedlist;
};

bool ValueCmpForErrorRecord(error_record const & a, error_record const & b)
{
    return a.cycle_number < b.cycle_number;
}

bool ValueCmpForInt(int const & a, int const & b)
{
    return a > b;
}
bool ValueCmpForLongInt(long int const & a, long int const & b)
{
    return a > b;
}

benchmark_listnode::benchmark_listnode(int s,int c,string n,int ran_round ,double *error_rate,int min, int max) {	
	signals = s;
	cycles = c;
	name = n;
	ran_times = ran_round;
	num_error = new int[ran_times];
	sqrt_win_size_min = min;
	sqrt_win_size_max = max;
	long int error_position_range = signals * cycles;
	//long int total_num_errors = error_position_range * (error_rate[0]/100);
	long int total_num_errors = cycles * (error_rate[0]/100);
	//cout << "total_num_errors: " << total_num_errors << endl;
	long int error_position_array[total_num_errors];

	error_record **array = new error_record*[ran_times];
	int j,k,l,rand_temp,rep,error_types_rand,num_signal_errors_rand_even,num_signal_errors_rand_odd,num_signal_errors;
	bool error_types_temp;
	int signal_number[signals];
	
	//result_record initialize
	int num_window_size = (log2(sqrt_win_size_max) - log2(sqrt_win_size_min) + 1);
	result_record_list_for_each_iteration = new result_record_linkedlist[NUMBER_OF_DIMENSIONS * num_window_size];

	for(j=0;j<num_window_size;j++) {
		for(k=0;k<NUMBER_OF_DIMENSIONS;k++) {
			int win_size_temp = pow(2,j*2) * sqrt_win_size_min * sqrt_win_size_min;
			int index = k + (j*NUMBER_OF_DIMENSIONS);
			result_record_list_for_each_iteration[index].SetWinSize(win_size_temp);
			result_record_list_for_each_iteration[index].SetDimension(k+1);
		}
	}
	
	//generating error
#ifdef DEBUG_EN
	cout << name << endl << "Signals: " << signals << endl << "Num_error: " << num_error << endl;;
#endif
	for(int i=0;i<ran_times;i++) {
#ifdef FAULT_INJECTION_FROM_FILE
		string circuit_name_fault_injected = "./signal_trace_data/signal_trace_data_" + n + "_FI_" + to_string(i) + "_list.txt";
		ifstream infile(circuit_name_fault_injected.c_str());
		int temp;
		string StringTemp;
		infile >> StringTemp;
		temp = stoi(StringTemp);
		num_error[i] = temp;
		array[i] = new error_record[num_error[i]];
		for(j=0;j<num_error[i];j++) {
			infile >> StringTemp;
			temp = stoi(StringTemp);
			array[i][j].cycle_number = temp;
			infile >> StringTemp;
			temp = stoi(StringTemp);
			num_signal_errors = temp;
			array[i][j].list.Set_num_error(num_signal_errors);
			for(l=0;l<num_signal_errors;l++) {
				infile >> StringTemp;
				temp = stoi(StringTemp);
				array[i][j].list.Push_front(temp);
			}
		}
#else
		num_error[i] = ceil(c * (error_rate[0]/100));
		array[i] = new error_record[num_error[i]];
		j=0;
		while(j < num_error[i]) {
			rand_temp = rand() % c;
			rep = 0;
			//Which cycles injecting fault
			for(k=0;k<j;k++) {
				if(rand_temp==array[i][k].cycle_number) {
					rep = 1;
				}			
			}
			//Which signals injecting fault
			if(rep==0) {
				if(signals==1)
					error_types_rand = rand() % (int) error_rate[1];
				else if(signals<=2)
					error_types_rand = rand() % (int) (error_rate[1]+error_rate[2]);
				else if(signals<=3)
					error_types_rand = rand() % (int) (error_rate[1]+error_rate[2]+error_rate[3]);
				else
					error_types_rand = rand() % 100;
				num_signal_errors_rand_odd = (signals%2) ? (rand() % (signals/2)) : (rand() % ((signals/2)-1));
				num_signal_errors_rand_even = rand() % ((signals/2)-1);

				if(error_types_rand < error_rate[1])
					num_signal_errors = 1;
				else if(error_types_rand < (error_rate[1]+error_rate[2]))
					num_signal_errors = 2;
				else if(error_types_rand < (error_rate[1]+error_rate[2]+error_rate[3]))
					num_signal_errors = 1 + ((num_signal_errors_rand_odd + 1)*2);
				else
					num_signal_errors = (num_signal_errors_rand_even + 2)*2;
				array[i][j].cycle_number = rand_temp;
				array[i][j].list.Set_num_error(num_signal_errors);
				//Shuffle
				for(k=0;k<signals;k++)
					signal_number[k] = k;
				for(k=signals-1;k>0;k--) {
					int rand_temp_l = rand()%k;
					int temp = signal_number[k];
					signal_number[k] = signal_number[rand_temp_l];
					signal_number[rand_temp_l] = temp;
				}
				sort(signal_number,signal_number+num_signal_errors, ValueCmpForInt);
				for(l=0;l<num_signal_errors;l++) {
					array[i][j].list.Push_front(signal_number[l]);
				}
				j++;
			}
		}
		sort(&array[i][0], &array[i][0]+num_error[i], ValueCmpForErrorRecord);

		//Random error generating
		/*j=0;
		while(j<total_num_errors) {
			long int error_position_rand = rand() % error_position_range;
			int re_rand_flag=0;
			for(l=0;l<j;l++) {
				if(error_position_rand==error_position_array[l]) {
					re_rand_flag = 1;
				}
			}
			if(re_rand_flag==0) {
				error_position_array[j] = error_position_rand;
				j++;
			}
		}
	
		sort(error_position_array,error_position_array+total_num_errors, ValueCmpForLongInt);
	
		long int num_cycle_errors=0;
		long int error_cycle_number=-1;
		//Computing number of cycle error
		for(j=0;j<total_num_errors;j++) {
			long int error_cycle_number_temp = error_position_array[j] / signals;
			if(error_cycle_number!=error_cycle_number_temp) {
				error_cycle_number = error_cycle_number_temp;
				num_cycle_errors++;
			}
		}
		num_error[i] = num_cycle_errors;
		array[i] = new error_record[num_error[i]];

		//Filled error_position_array iteration init
		long int num_signal_errors=1;
		long int error_cycle_index = 0;
		int error_signal_number=0;
	
		//Filled error_position_array[0]
		error_signal_number = error_position_array[0] % signals;
		error_cycle_number = error_position_array[0] / signals;
		array[i][error_cycle_index].list.Push_front(error_signal_number);

		for (j=1;j<total_num_errors;j++) {
			long int error_cycle_number_temp = error_position_array[j] / signals;
			error_signal_number = error_position_array[j] % signals;
			if(error_cycle_number==error_cycle_number_temp) {
				num_signal_errors++;
				array[i][error_cycle_index].list.Push_front(error_signal_number);
			}
			else {
				array[i][error_cycle_index].list.Set_num_error(num_signal_errors);
				array[i][error_cycle_index].cycle_number = error_cycle_number;
				error_cycle_index++;
				array[i][error_cycle_index].list.Push_front(error_signal_number);

				error_cycle_number = error_cycle_number_temp;
				num_signal_errors = 1;
			}
		
		}
		array[i][error_cycle_index].list.Set_num_error(num_signal_errors);
		array[i][error_cycle_index].cycle_number = error_cycle_number;

		sort(&array[i][0], &array[i][0]+num_cycle_errors, ValueCmpForErrorRecord);*/
#endif

#ifdef WRITE_FAULT_INJECTED_RESULT
		string circuit_name_fault_injected = "./signal_trace_data/signal_trace_data_" + n + "_FI_" + to_string(i) + "_list.txt";
		ofstream outfile(circuit_name_fault_injected.c_str());

		outfile << num_error[i] << endl;
		for(j=0;j<num_error[i];j++) {
			outfile << array[i][j].cycle_number << endl;
			int num_signal_errors_temp = array[i][j].list.Get_num_error();
			outfile << num_signal_errors_temp << endl;
			for(k=(num_signal_errors_temp-1);k>=0;k--)
				outfile << array[i][j].list.Get_signal_number(k) << " ";
			outfile << endl;
		}
#endif
	}
	errors = array;
}

void benchmark_listnode::DisplayErrors(){
	int i,j,k,num_signal_error;
	error_record_for_signals temp;
	//cout << name << "\tnumber of error: " << num_error << endl;
	for (i=0;i<ran_times;i++) {
		cout << "ran_" << i << endl;	
		for(j=0;j<num_error[i];j++) {
			num_signal_error = errors[i][j].list.Get_num_error();
			cout << "Cycle number: " << errors[i][j].cycle_number << "\tnumber of signal error " << num_signal_error << endl;
			for(k=0;k<num_signal_error;k++)	{
				cout << errors[i][j].list.Get_signal_number(k) << "\t";
				if((k%5)==4) cout << endl;	
			}
			cout << endl;
		}
		cout << endl;
	}
}

void benchmark_listnode::DisplayErrors(int ran_time_number){
	int i,j,k,num_signal_error;
	error_record_for_signals temp;
	cout << name << "\tnumber of error: " << num_error[ran_time_number] << endl;	
	for(j=0;j<num_error[ran_time_number];j++) {
		num_signal_error = errors[ran_time_number][j].list.Get_num_error();
		cout << "Cycle number: " << errors[ran_time_number][j].cycle_number << "\tnumber of signal error " << num_signal_error << endl;
		for(k=0;k<num_signal_error;k++)	{
			cout << errors[ran_time_number][j].list.Get_signal_number(k) << "\t";
			if((k%5)==4) cout << endl;	
		}
		cout << endl;
	}
	cout << endl;
}

int benchmark_listnode::GetNumError(int ran_time_number){
	return num_error[ran_time_number];
}

int benchmark_listnode::GetSignals(){
	return signals;
}

int benchmark_listnode::GetCycles(){
	return cycles;
}

int benchmark_listnode::GetRanTimes(){
	return ran_times;
}

string benchmark_listnode::GetName(){
	return name;
}

int benchmark_listnode::GetWinSizeMin(){
	return sqrt_win_size_min;
}

int benchmark_listnode::GetWinSizeMax(){
	return sqrt_win_size_max;
}


class benchmark_linkedlist {
	private:
		benchmark_listnode *first;
	public:
		benchmark_linkedlist():first(0) {};
		void PrintList();
		void PrintErrors();
    		void Push_front(int s,int c,string n,int ran_round ,double *error_rate,int min, int max);
		benchmark_listnode Pop_front();
		benchmark_listnode Get_benchmark_listnode(int index);
    		void Delete(string n);
    		void Clear();
};

void benchmark_linkedlist::PrintList(){

    if (first == NULL) {                      
        cout << "List is empty.\n";
        return;
    }

    benchmark_listnode *current = first;            
    while (current != NULL) {                 
        cout << current->name << "\t" << current->signals << "\t" << current->cycles << endl;
        current = current->next;
    }
    cout << endl;
}

void benchmark_linkedlist::PrintErrors(){

    if (first == NULL) {                      
        cout << "List is empty.\n";
        return;
    }

    benchmark_listnode *current = first;            
    while (current != NULL) {                 
        current->DisplayErrors();
        current = current->next;
    }
    cout << endl;
}

void benchmark_linkedlist::Push_front(int s,int c,string n,int ran_round ,double *error_rate,int min, int max){

    benchmark_listnode *newNode = new benchmark_listnode(s,c,n,ran_round,error_rate,min,max);   
    newNode->next = first;                 
    first = newNode;                       
}

benchmark_listnode benchmark_linkedlist::Pop_front(){

    benchmark_listnode IndexNode = *first;
    first = first->next;
    return IndexNode;                      
}

benchmark_listnode benchmark_linkedlist::Get_benchmark_listnode(int index){
    int i;
    benchmark_listnode popNode = *first;
    for(i=0;i<index;i++)
	popNode = *(popNode.next);
    return popNode;                      
}

void benchmark_linkedlist::Delete(string n){

    benchmark_listnode *current = first, *previous = 0;

    while (current != 0 && current->name != n) {  
        previous = current;                       
        current = current->next;                  
    }                                            

    if (current == 0) {                
        cout << "There is no " << n << " in list.\n";
    }
    else if (current == first) {        
        first = current->next;
        delete current;
        current = 0;
    }
    else {                               
        previous->next = current->next; 
        delete current;
        current = 0;
    }
}

void benchmark_linkedlist::Clear(){
    while (first != 0) {            // Traversal
        benchmark_listnode *current = first;
        first = first->next;
	delete current->errors;
        delete current;
        current = 0;
    }
}
