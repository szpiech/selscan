# Makefile for selscan

#G++FLAG = -g -Wall -Werror   #debugging

#for osx/macos-arm systems
CC = g++-14 -g -std=c++17
#CC = g++-13 -g
G++FLAG = -O3 -m64  #-mmmx -msse -msse2
LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas #For dynamic linking
#LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a   #For static linking of norm program to gsl libs
LINK_OPTS2 = -lpthread -lz -fopenmp 
I_PATH = -I../include
L_PATH = ../lib/macos-arm #L_PATH = ../lib/osx

OBJ_STATS = xpihh.o ehh.o ihs.o

all : selscan norm
# build_and_cleanup

# build_and_cleanup: selscan norm
# 	rm -rf *.o	

selscan : selscan-main.o selscan-maintools.o $(OBJ_STATS) selscan-data.o selscan-cli.o binom.o selscan-pbar.o param_t.o gzstream.o hamming_t.o  
	$(CC) $(G++FLAG) $(LINK_OPTS2) -o selscan selscan-main.o selscan-maintools.o $(OBJ_STATS) selscan-data.o selscan-cli.o binom.o selscan-pbar.o param_t.o gzstream.o hamming_t.o
	
selscan-main.o : selscan-main.cpp
	$(CC) $(G++FLAG) -c selscan-main.cpp $(I_PATH)

selscan-maintools.o : selscan-maintools.cpp selscan-maintools.h 
	$(CC) $(G++FLAG) -c selscan-maintools.cpp $(I_PATH)

%.o: stats/%.cpp stats/%.h 
	$(CC) $(G++FLAG) -c $< -o $@ $(I_PATH)

selscan-data.o : selscan-data.cpp selscan-data.h
	$(CC) $(G++FLAG) -c selscan-data.cpp $(I_PATH)

selscan-cli.o : selscan-cli.cpp selscan-cli.h
	$(CC) $(G++FLAG) -c selscan-cli.cpp $(I_PATH)

selscan-pbar.o : selscan-pbar.cpp selscan-pbar.h
	$(CC) $(G++FLAG) -c selscan-pbar.cpp $(I_PATH)

binom.o : binom.cpp binom.h
	$(CC) $(G++FLAG) -c binom.cpp

norm : norm.o param_t.o
	$(CC) $(G++FLAG) -o norm norm.o param_t.o $(LINK_OPTS)

norm.o : norm.cpp
	$(CC) $(G++FLAG) -c norm.cpp $(I_PATH)

param_t.o : param_t.cpp param_t.h
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp gzstream.h
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)

hamming_t.o : hamming_t.cpp hamming_t.h
	$(CC) $(G++FLAG) -c hamming_t.cpp $(I_PATH)

# $(BIN_DIR): # Create the binary directory if it does not exist
# 	mkdir -p $(BIN_DIR)

.PHONY: cleanup
.PHONY: clean

clean :
	rm -rf *.o;  rm -rf selscan; rm -rf  norm; rm -rf outfile*
