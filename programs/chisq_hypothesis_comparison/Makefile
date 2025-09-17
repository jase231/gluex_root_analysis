CXX = g++
CXXFLAGS = -O -Wall
ROOTFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS := $(shell $(ROOTSYS)/bin/root-config --libs)

# Source files
SRCS = compare_hypotheses.cpp main.cpp
# Object files
OBJS = $(SRCS:.cpp=.o)
# Executable name
EXEC = compare_hypotheses

# Default rule: build the executable (release mode)
all: $(EXEC)

# Rule to link the executable
$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Rule to compile source files into object files
%.o: %.cpp
	$(CXX) -c $(ROOTFLAGS) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean
