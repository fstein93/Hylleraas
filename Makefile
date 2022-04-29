# Thanks to Job Vranish (https://spin.atomicobject.com/2016/08/26/makefile-c-projects/)
_TARGET_EXEC := Hylleraas

BUILD_DIR := ./build/release
SRC_DIR := ./src
EXE_RELEASE_DIR := ./exe/release
EXE_DEBUG_DIR := ./exe/debug

LDFLAGS := -llapack -lblas -lm
INCFLAGS := -I/usr/include/x86_64-linux-gnu/cblas.h
LIBS := -L/usr/lib/x86_64-linux-gnu/ -L/usr/lib/x86_64-linux-gnu/blas

EXE_RELEASE := $(EXE_RELEASE_DIR)/$(_TARGET_EXEC)
EXE_DEBUG := $(EXE_DEBUG_DIR)/$(_TARGET_EXEC)

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. Make will incorrectly expand these otherwise.
SRCS := $(shell find $(SRC_DIR) -name '*.cpp')

# String substitution for every C/C++ file.
# As an example, hello.cpp turns into ./build/hello.cpp.o
OBJS_RELEASE := $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIR) -type d) $(LIBS) $(INCFLAGS)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := $(INC_FLAGS) -MMD -MP

all: $(EXE_RELEASE) $(EXE_DEBUG)

# The final build step.
$(EXE_RELEASE): $(OBJS_RELEASE)
	mkdir -p $(dir $@)
	$(CXX) $(OBJS_RELEASE) -o $@ $(LDFLAGS)

$(EXE_DEBUG): $(OBJS_RELEASE)
	mkdir -p $(dir $@)
	$(CXX) $(OBJS_RELEASE) -o $@ $(LDFLAGS)

# Build step for C++ source
$(OBJS_RELEASE): $(SRCS)
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(EXE_RELEASE_DIR) $(EXE_DEBUG_DIR)

