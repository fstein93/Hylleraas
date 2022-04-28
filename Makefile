# Project: Projekt1
# Makefile created by Dev-C++ 5.11

CPP      = g++
CC       = gcc

_SRC = Hylleraas.cpp
_DEPS = integrator.h


SRC = $(_SRC)
DEPS = $(_DEPS)

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
	$(CPP) -o $@ $< $(LIBS) $(RELEASEFLAGS)

$(BIN_DEBUG): $(OBJ_DEBUG)
	$(CPP) -o $@ $< $(LIBS) $(DEBUGFLAGS)

$(OBJ): $(SRC) $(DEPS)
	$(CPP) -c -o $@ $< $(CXXFLAGS) $(RELEASEFLAGS)

$(OBJ_DEBUG): $(SRC) $(DEPS)
	$(CPP) -c -o $@ $< $(CXXFLAGS) $(DEBUGFLAGS)

$(DEPS):
	

