EVOSIM_ROOT = $(shell pwd)
export PATH := $(shell pwd):$(PATH)
BINDIR = $(EVOSIM_ROOT)/bin

ifndef SMITHLAB_CPP
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/src/smithlab_cpp
ifeq ("$(wildcard $(SMITHLAB_CPP))","")
$(error SMITHLAB_CPP variable not set and smithlab_cpp not found)
endif
endif

ifndef TREETOOL
TREETOOL=$(abspath $(dir $(MAKEFILE_LIST)))/src/adssrc/treetool
ifeq ("$(wildcard $(TREETOOL))","")
$(error could not set TREETOOL variable)
endif
endif


all:
	@make -C src EVOSIM_ROOT=$(EVOSIM_ROOT) SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) OPT=1

install:
	@export PATH=$(PATH)
	@make -C src EVOSIM_ROOT=$(EVOSIM_ROOT) SMITHLAB_CPP=$(SMITHLAB_CPP) TREETOOL=$(TREETOOL) OPT=1 install

test:
	@make -C src OPT=1 test
.PHONY: test 

clean:
	@rm -rf $(BINDIR)
	@make -C src EVOSIM_ROOT=$(EVOSIM_ROOT) clean
.PHONY: clean
