# Makefile for selscan

#G++FLAG = -g -Wall -Werror
G++FLAG = -O3 -m64 -mmmx -msse -msse2
I_PATH = ../include
L_PATH = ../lib

#For static linking
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a

#For dynamic linking
#LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

all : selscan norm

selscan : selscan-main.o selscan-data.o binom.o selscan-pbar.o param_t.o
	g++ -o selscan selscan-main.o selscan-data.o binom.o selscan-pbar.o param_t.o -pthread

selscan-main.o : selscan-main.cpp
	g++ $(G++FLAG) -c selscan-main.cpp 

selscan-data.o : selscan-data.cpp
	g++ $(G++FLAG) -c selscan-data.cpp 

selscan-pbar.o : selscan-pbar.cpp
	g++ $(G++FLAG) -c selscan-pbar.cpp 

binom.o : binom.cpp
	g++ $(G++FLAG) -c binom.cpp

norm : norm.o param_t.o
	g++ -o norm norm.o param_t.o $(LINK_OPTS)

norm.o : norm.cpp
	g++ $(G++FLAG) -c norm.cpp -I$(I_PATH)

param_t.o : param_t.cpp
	g++ $(G++FLAG) -c param_t.cpp

clean :
	rm *.o
