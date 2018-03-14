# Copyright (C) 2018 University of Southern California
#                    Jianghan Qu and Andrew D Smith
#
# Author: Andrew D. Smith and Jianghan Qu
#
# This is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA

PROGS = evosim estparam epievo_sim global_jumps_to_paths \
	extract_state_sequences test_sampling

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50 -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g
LIBS = -lgsl -lgslcblas

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

COMMON_DIR = $(realpath ../common)
INCLUDEDIRS = $(COMMON_DIR) $(SMITHLAB_CPP) $(TREETOOL)
INCLUDEARGS = $(addprefix -I, $(INCLUDEDIRS))

all: $(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, OptionParser.o smithlab_utils.o smithlab_os.o) \
	$(addprefix $(COMMON_DIR)/, PhyloTreePreorder.o TripletPattern.o param.o jump.o Path.o) \
	$(addprefix $(TREETOOL)/, PhyloTree.o )

epievo_sim: $(addprefix $(COMMON_DIR)/, TripletSampler.o StateSeq.o \
	EpiEvoModel.o GlobalJump.o)

test_sampling: $(addprefix $(COMMON_DIR)/, TripletSampler.o StateSeq.o \
	EpiEvoModel.o)

global_jumps_to_paths: $(addprefix $(COMMON_DIR)/, GlobalJump.o StateSeq.o)

extract_state_sequences: $(addprefix $(COMMON_DIR)/, TripletSampler.o \
	StateSeq.o GlobalJump.o)


install: all
	@mkdir -p ../../bin
	@install -m 755 $(PROGS) ../../bin

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
