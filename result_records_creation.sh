#/bin/tcsh
set bug_rate_list="0.01 0.05 0.1"
set exploration_template_file="./exploration_tool_template.cpp"
set exploration_file="./exploration_tool.cpp"
set exploration_exe="exploration_tool"
set result_records_file="./Result_records/Result_records_voting.txt"
set result_records_lib_dir="./Result_records_lib_PG_N4/"

mkdir -p $result_records_lib_dir

foreach brl ($bug_rate_list)
	#Inject SEU bugs
	cp $exploration_template_file $exploration_file
	sed -i "s/error_rate\[1\] = N/error_rate[1] = 100/g" $exploration_file
	sed -i "s/error_rate\[2\] = N/error_rate[2] = 0/g" $exploration_file
	sed -i "s/error_rate\[3\] = N/error_rate[3] = 0/g" $exploration_file
	sed -i "s/error_rate\[4\] = N/error_rate[4] = 0/g" $exploration_file
	sed -i "s/error_rate\[0\] = N/error_rate[0] = $brl/g" $exploration_file
	make
	time ./$exploration_exe > t.log
	set result_records_lib_file="$brl""_SEU.txt"
	mv $result_records_file "$result_records_lib_dir$result_records_lib_file"
	#break
	#Inject Random bugs
	cp $exploration_template_file $exploration_file
	sed -i "s/error_rate\[1\] = N/error_rate[1] = 25/g" $exploration_file
	sed -i "s/error_rate\[2\] = N/error_rate[2] = 25/g" $exploration_file
	sed -i "s/error_rate\[3\] = N/error_rate[3] = 25/g" $exploration_file
	sed -i "s/error_rate\[4\] = N/error_rate[4] = 25/g" $exploration_file
	sed -i "s/error_rate\[0\] = N/error_rate[0] = $brl/g" $exploration_file
	make
	time ./$exploration_exe > t.log
	set result_records_lib_file="$brl""_RAN.txt"
	mv $result_records_file "$result_records_lib_dir$result_records_lib_file"
end
