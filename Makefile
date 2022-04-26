# Project: Projekt1
# Makefile created by Dev-C++ 5.11

CPP      = g++
CC       = gcc
WINDRES  = windres.exe
OBJ      = Hylleraas.o
OBJ_DEBUG = Hylleraas.o_debug
LIBS     = 
INCS     = 
CXXINCS  = 
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

Hylleraas.o: Hylleraas.cpp
	$(CPP) -c Hylleraas.cpp -o Hylleraas.o $(CXXFLAGS) $(RELEASEFLAGS)

Hylleraas.o_debug: Hylleraas.cpp
	$(CPP) -c Hylleraas.cpp -o Hylleraas.o_debug $(CXXFLAGS) $(DEBUGFLAGS)
