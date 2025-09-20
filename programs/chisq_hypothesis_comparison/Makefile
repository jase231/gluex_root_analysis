CXX = g++
CXXFLAGS = -O -Wall
ROOTFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS = $(shell $(ROOTSYS)/bin/root-config --libs)

SRCS = compare_hypotheses.cpp main.cpp
EXEC = chisq_hypothesis_comparison
BIN_DIR = $(BMS_OSNAME)/bin
BIN_PATH = $(BIN_DIR)/$(EXEC)
OBJS = $(SRCS:%.cpp=$(BIN_DIR)/%.o)

all: install

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

# Rule to compile source files into object files
$(BIN_DIR)/%.o: %.cpp | $(BIN_DIR)
	$(CXX) -c $(ROOTFLAGS) $(CXXFLAGS) $< -o $@

# Rule to link the executable
$(BIN_PATH): $(OBJS) | $(BIN_DIR)
	@rm -f $(BIN_PATH)
	$(CXX) -o $(BIN_PATH) $(OBJS) $(LDFLAGS)
	@chmod 755 $(BIN_PATH)

install: $(BIN_PATH)
	@mkdir -p ../../$(BIN_DIR)
	@cp $(BIN_PATH) ../../$(BIN_DIR)/

clean:
	rm -rf $(BMS_OSNAME)
	rm -f ../../$(BIN_PATH)

env:
	@echo BMS_OSNAME $(BMS_OSNAME)
	@echo FLAGS $(CXXFLAGS) $(ROOTFLAGS)
	@echo LIBS $(LDFLAGS)

.PHONY: all clean install env
