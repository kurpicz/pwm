################################################################################
# Makefile
#
# Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

CFLAGS := -ansi
CXXFLAGS := -std=c++14 -O3 -ffast-math -funroll-loops -march=native -fopenmp

LIBS := -ldl

SRC_DIR := benchmark/src
BIN_DIR := benchmark/bin
INC_DIR := include

MC_DIR := external/malloc_count
MC_OBJS := $(MC_DIR)/malloc_count.o

BENCHMARKS := construct check memory

WM_ALGORITHS := $(foreach ALGORITHM, $(wildcard $(INC_DIR)/wm_*.hpp),\
	$(foreach BENCHMARK, $(BENCHMARKS),\
		$(BIN_DIR)/wm/$(BENCHMARK)_$(basename $(notdir $(ALGORITHM)))))

WT_ALGORITHS := $(foreach ALGORITHM, $(wildcard $(INC_DIR)/wt_*.hpp),\
	$(foreach BENCHMARK, $(BENCHMARKS),\
		$(BIN_DIR)/wt/$(BENCHMARK)_$(basename $(notdir $(ALGORITHM)))))

execs: $(WM_ALGORITHS) $(WT_ALGORITHS)

$(BIN_DIR)/wm/construct_%: $(SRC_DIR)/construct_wm.cpp $(INC_DIR)/%.hpp
	@mkdir -p $(BIN_DIR)/wm
	@echo "Compiling construct_$*"
	@$(CXX) $(CXXFLAGS) -DWM_TYPE="$*" -DRUNS=5 -DTIMING=1\
		$(SRC_DIR)/construct_wm.cpp -I${INC_DIR} -o $@

$(BIN_DIR)/wm/check_%: $(SRC_DIR)/construct_wm.cpp $(INC_DIR)/%.hpp
	@mkdir -p $(BIN_DIR)/wm
	@echo "Compiling check_$*"
	@$(CXX) $(CXXFLAGS) -DWM_TYPE="$*" -DRUNS=1 -DCHECK=1\
		$(SRC_DIR)/construct_wm.cpp -I${INC_DIR} -o $@

$(BIN_DIR)/wm/memory_%: $(SRC_DIR)/construct_wm.cpp $(INC_DIR)/%.hpp $(MC_OBJS)
	@mkdir -p $(BIN_DIR)/wm
	@echo "Compiling check_$*"
	@$(CXX) $(CXXFLAGS) -DWM_TYPE="$*" -DRUNS=1 -DMEMORY=1\
		$(SRC_DIR)/construct_wm.cpp -I${INC_DIR} -o $@ $(MC_OBJS) $(LIBS)

$(BIN_DIR)/wt/construct_%: $(SRC_DIR)/construct_wt.cpp $(INC_DIR)/%.hpp
	@mkdir -p $(BIN_DIR)/wt
	@echo "Compiling construct_$*"
	@$(CXX) $(CXXFLAGS) -DWT_TYPE="$*" -DRUNS=5 -DTIMING=1\
		$(SRC_DIR)/construct_wt.cpp -I${INC_DIR} -o $@

$(BIN_DIR)/wt/check_%: $(SRC_DIR)/construct_wt.cpp $(INC_DIR)/%.hpp
	@mkdir -p $(BIN_DIR)/wt
	@echo "Compiling check_$*"
	@$(CXX) $(CXXFLAGS) -DWT_TYPE="$*" -DRUNS=1 -DCHECK=1\
		$(SRC_DIR)/construct_wt.cpp -I${INC_DIR} -o $@

$(BIN_DIR)/wt/memory_%: $(SRC_DIR)/construct_wt.cpp $(INC_DIR)/%.hpp $(MC_OBJS)
	@mkdir -p $(BIN_DIR)/wt
	@echo "Compiling check_$*"
	@$(CXX) $(CXXFLAGS) -DWT_TYPE="$*" -DRUNS=1 -DMEMORY=1\
		$(SRC_DIR)/construct_wt.cpp -I${INC_DIR} -o $@ $(MC_OBJS) $(LIBS)

$(MC_DIR)/%.o : $(MC_DIR)/%c
	$(CC) $(CFLAGS) -o $@ $<

.PHONY: clean
clean:
	@rm -f $(BIN_DIR)/wm/*
	@rm -f $(BIN_DIR)/wt/*
