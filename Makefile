# Project: Projekt1
# Makefile created by Dev-C++ 5.11

CPP      = g++
CC       = gcc
OBJ      = Hylleraas.o
OBJ_DEBUG = Hylleraas.o_debug
LIBS     = -llapack -lblas
INCS     = -I/usr/include/x86_64-linux-gnu/cblas.h
CXXINCS  = -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/blas
BIN      = Hylleraas.exe
BIN_DEBUG = Hylleraas_debug.exe
CXXFLAGS = $(CXXINCS) -march=native -std=gnu++11 -Wall -Wextra -pedantic
CFLAGS   = $(INCS) -march=native -std=gnu++11 -Wall -Wextra -pedantic
DEBUGFLAGS = -pg -Og -g
RELEASEFLAGS = -O2
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) $(BIN_DEBUG) all-after

clean: clean-custom
	${RM} $(OBJ) $(OBJ_DEBUG) $(BIN) $(BIN_DEBUG)

$(BIN): $(OBJ)
	$(CPP) $(OBJ) -o $(BIN) $(LIBS) $(RELEASEFLAGS)

$(BIN_DEBUG): $(OBJ_DEBUG)
	$(CPP) $(OBJ_DEBUG) -o $(BIN_DEBUG) $(LIBS) $(DEBUGFLAGS)

$(OBJ): Hylleraas.cpp
	$(CPP) -c Hylleraas.cpp -o Hylleraas.o $(CXXFLAGS) $(RELEASEFLAGS)

$(OBJ_DEBUG): Hylleraas.cpp
	$(CPP) -c Hylleraas.cpp -o Hylleraas.o_debug $(CXXFLAGS) $(DEBUGFLAGS)
