PROG=main

PROG_OBJ= binary_bcast_one_side.o binary_bcast.o binomial_bcast_one_sided.o binomial_bcast.o linear_bcast.o $(PROG).o
INCLUDE_PATH=./

CXX=mpic++
CXXFLAGS=-Wall -std=c++17 -I$(INCLUDE_PATH) -pthread
#		 -lboost_system -lboost_date_time -lboost_thread \
#		 -L/usr/lib64

all: $(PROG)

$(PROG): $(PROG_OBJ) 
	$(CXX) $(CXXFLAGS) $^ -o $@ 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

clean:
	rm -rf $(PROG) $(PROG_OBJ) 
	rm -rf stdout
	rm -rf mpitask.o*
	rm -rf *.out