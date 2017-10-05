CC = g++ -std=c++0x
OBJ = exploration_tool.o
APP = exploration_tool
OUT_DIR = DecisionModeOutput Result_records

all: directories exploration_tool

directories: $(OUT_DIR)

${OUT_DIR}:
	mkdir -p $(OUT_DIR)	
	
%.o: %.cpp
	$(CC) -c -o $@ $< -pthread

exploration_tool: $(OBJ)
	$(CC) -o $@ $^ -pthread

clear:
	rm -f *.o $(APP) 
	rm -rf $(OUT_DIR)
	rm -f ./signal_trace_data/*_FI_*
