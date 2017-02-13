################################################################################
# Makefile
#
# Copyright (C) 2017 Florian Kurpicz <florian.kurpicz@tu-dortmund.de>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################

CFLAGS := -std=c++14 -O3 -march=native -DNDEBUG

SRC_DIR := benchmark/src
BIN_DIR := benchmark/bin
INC_DIR := include

BENCHMARKS := construct check

WM_ALGORITHS = $(foreach ALGORITHM, $(wildcard $(INC_DIR)/wm_*.hpp),\
	$(foreach BENCHMARK, $(BENCHMARKS),\
		$(BIN_DIR)/$(BENCHMARK)_$(basename $(notdir $(ALGORITHM)))))

execs: $(WM_ALGORITHS)

$(BIN_DIR)/construct_%: $(SRC_DIR)/construct_wm.cpp $(INC_DIR)/%.hpp
	@echo "Compiling construct_$*"
	$(CXX) $(CFLAGS) -DWM_TYPE="$*" -DRUNS=5 -DTIMING=1\
		$(SRC_DIR)/construct_wm.cpp -I${INC_DIR} -o $@

$(BIN_DIR)/check_%: $(SRC_DIR)/construct_wm.cpp $(INC_DIR)/%.hpp
	@echo "Compiling check_$*"
	$(CXX) $(CFLAGS) -DWM_TYPE="$*" -DRUNS=1 -DCHECK=1\
		$(SRC_DIR)/construct_wm.cpp -I${INC_DIR} -o $@

.PHONY: clean
clean:
	@rm -f $(BIN_DIR)/*