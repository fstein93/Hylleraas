# Thanks to Job Vranish (https://spin.atomicobject.com/2016/08/26/makefile-c-projects/)
_TARGET_EXEC := Hylleraas

BUILD_RELEASE_DIR := ./build/release
BUILD_DEBUG_DIR := ./build/debug
SRC_DIR := ./src
EXE_RELEASE_DIR := ./exe/release
EXE_DEBUG_DIR := ./exe/debug

LDFLAGS := -llapack -lblas -lm
INCFLAGS := -I/usr/include/x86_64-linux-gnu/
LIBS := -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/blas
WARN_FLAGS :=  -Wall -Werror -Wextra -pedantic -Wshadow -Wsign-conversion -Wunreachable-code -Wconversion
RELEASE_FLAGS := -O2 -g0 -march=native -std=c++11 -fdelete-dead-exceptions
DEBUG_FLAGS := -Og -g -pg -march=native -std=c++11 -fbounds-check -fstack-protector -fcf-protection -fsanitize=leak,undefined
# -fsanitize=address

EXE_RELEASE := $(EXE_RELEASE_DIR)/$(_TARGET_EXEC)
EXE_DEBUG := $(EXE_DEBUG_DIR)/$(_TARGET_EXEC)

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. Make will incorrectly expand these otherwise.
SRCS := $(shell find $(SRC_DIR) -name '*.cpp')
HEADERS := $(shell find $(SRC_DIR) -name '*.h')

# String substitution for every C/C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS_RELEASE := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_RELEASE_DIR)/%.o)
OBJS_DEBUG := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DEBUG_DIR)/%.o)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(LIBS) $(INCFLAGS)
INC_FLAGS := $(INC_DIRS)

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := $(INC_FLAGS) $(WARN_FLAGS) -MMD -MP

all: $(EXE_RELEASE) $(EXE_DEBUG)

# The final build step.
$(EXE_RELEASE): $(OBJS_RELEASE)
	mkdir -p $(dir $@)
	$(CXX) $(OBJS_RELEASE) -o $@ $(LDFLAGS) $(RELEASE_FLAGS) $(WARN_FLAGS)

$(EXE_DEBUG): $(OBJS_DEBUG)
	mkdir -p $(dir $@)
	$(CXX) $(OBJS_DEBUG) -o $@ $(LDFLAGS) $(DEBUG_FLAGS) $(WARN_FLAGS)

# Build step for C++ source
$(OBJS_RELEASE): $(SRCS) $(HEADERS)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(RELEASE_FLAGS) -c $< -o $@

# Build step for C++ source
$(OBJS_DEBUG): $(SRCS) $(HEADERS)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEBUG_FLAGS) -c $< -o $@

check: $(EXE_DEBUG)
	$(EXE_DEBUG) < test.inp

.PHONY: clean
clean:
	rm -rf $(BUILD_RELEASE_DIR) $(BUILD_DEBUG_DIR) $(EXE_RELEASE_DIR) $(EXE_DEBUG_DIR)

