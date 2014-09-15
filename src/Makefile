# Makefile for selscan

#debugging
#G++FLAG = -g -Wall -Werror
#For dynamic linking
#LINK_OPTS = -L$(L_PATH) -lgsl -lgslcblas

#for osx systems
#CC = g++
#G++FLAG = -fast -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
#I_PATH = -I../include
#L_PATH = ../lib/osx

#for linux systems
#CC = g++
#G++FLAG = -O3 -m64 -mmmx -msse -msse2
#LINK_OPTS2 = -pthread -lz
#I_PATH = -I../include
#L_PATH = ../lib/linux

#for windows, using MinGW build environment
CC = g++.exe
G++FLAG = -DPTW32_STATIC_LIB -O3 -static-libgcc -static-libstdc++
LINK_OPTS2 = ../lib/win32/libpthreadGC2.a ../lib/win32/libz.a 
I_PATH = -I../include -I../include/win32 
L_PATH = ../lib/win32

#For static linking of norm program to gsl libs
LINK_OPTS = $(L_PATH)/libgsl.a $(L_PATH)/libgslcblas.a

all : selscan norm

selscan : selscan-main.o selscan-data.o binom.o selscan-pbar.o param_t.o gzstream.o
	$(CC) $(G++FLAG) -o selscan selscan-main.o selscan-data.o binom.o selscan-pbar.o param_t.o gzstream.o $(LINK_OPTS2)

selscan-main.o : selscan-main.cpp
	$(CC) $(G++FLAG) -c selscan-main.cpp $(I_PATH)

selscan-data.o : selscan-data.cpp
	$(CC) $(G++FLAG) -c selscan-data.cpp $(I_PATH)

selscan-pbar.o : selscan-pbar.cpp
	$(CC) $(G++FLAG) -c selscan-pbar.cpp $(I_PATH)

binom.o : binom.cpp
	$(CC) $(G++FLAG) -c binom.cpp

norm : norm.o param_t.o
	$(CC) $(G++FLAG) -o norm norm.o param_t.o $(LINK_OPTS)

norm.o : norm.cpp
	$(CC) $(G++FLAG) -c norm.cpp $(I_PATH)

param_t.o : param_t.cpp
	$(CC) $(G++FLAG) -c param_t.cpp

gzstream.o : gzstream.cpp
	$(CC) $(G++FLAG) -c gzstream.cpp $(I_PATH)

clean :
	rm *.o
