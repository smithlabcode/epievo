# This file is part of epievo
#
# Copyright (C) 2020: University of Southern California
#                     Andrew D. Smith and Xiaojing Ji
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I $(top_srcdir)/src/libepievo -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 # default has optimization on

noinst_LIBRARIES = libepievo.a libsmithlab_cpp.a
libepievo_a_SOURCES = \
	src/libepievo/ContinuousTimeMarkovModel.cpp \
	src/libepievo/EndCondSampling.cpp \
	src/libepievo/EpiEvoModel.cpp \
	src/libepievo/GlobalJump.cpp \
	src/libepievo/IndepSite.cpp \
	src/libepievo/IntervalSampler.cpp \
	src/libepievo/ParamEstimation.cpp \
	src/libepievo/Path.cpp \
	src/libepievo/PhyloTree.cpp \
	src/libepievo/PhyloTreePreorder.cpp \
	src/libepievo/Segment.cpp \
	src/libepievo/SingleSiteSampler.cpp \
	src/libepievo/TreeHelper.cpp \
	src/libepievo/TripletSampler.cpp \
	src/libepievo/epievo_utils.cpp

libepievo_a_SOURCES += \
	src/libepievo/ContinuousTimeMarkovModel.hpp \
	src/libepievo/EndCondSampling.hpp \
	src/libepievo/EpiEvoModel.hpp \
	src/libepievo/GlobalJump.hpp \
	src/libepievo/IndepSite.hpp \
	src/libepievo/IntervalSampler.hpp \
	src/libepievo/ParamEstimation.hpp \
	src/libepievo/Path.hpp \
	src/libepievo/PhyloTree.hpp \
	src/libepievo/PhyloTreePreorder.hpp \
	src/libepievo/Segment.hpp \
	src/libepievo/SingleSiteSampler.hpp \
	src/libepievo/TreeHelper.hpp \
	src/libepievo/TripletSampler.hpp \
	src/libepievo/emission_utils.hpp \
	src/libepievo/epievo_utils.hpp

libsmithlab_cpp_a_SOURCES = \
	src/smithlab_cpp/GenomicRegion.cpp \
	src/smithlab_cpp/MappedRead.cpp \
	src/smithlab_cpp/OptionParser.cpp \
	src/smithlab_cpp/QualityScore.cpp \
	src/smithlab_cpp/bisulfite_utils.cpp \
	src/smithlab_cpp/chromosome_utils.cpp \
	src/smithlab_cpp/dna_four_bit.cpp \
	src/smithlab_cpp/htslib_wrapper.cpp \
	src/smithlab_cpp/sim_utils.cpp \
	src/smithlab_cpp/smithlab_os.cpp \
	src/smithlab_cpp/smithlab_utils.cpp \
	src/smithlab_cpp/zlib_wrapper.cpp

libsmithlab_cpp_a_SOURCES += \
	src/smithlab_cpp/GenomicRegion.hpp \
	src/smithlab_cpp/MappedRead.hpp \
	src/smithlab_cpp/OptionParser.hpp \
	src/smithlab_cpp/QualityScore.hpp \
	src/smithlab_cpp/bisulfite_utils.hpp \
	src/smithlab_cpp/chromosome_utils.hpp \
	src/smithlab_cpp/cigar_utils.hpp \
	src/smithlab_cpp/dna_four_bit.hpp \
	src/smithlab_cpp/htslib_wrapper.hpp \
	src/smithlab_cpp/sim_utils.hpp \
	src/smithlab_cpp/smithlab_os.hpp \
	src/smithlab_cpp/smithlab_utils.hpp \
	src/smithlab_cpp/zlib_wrapper.hpp

LDADD = libepievo.a libsmithlab_cpp.a

### PROG SUBDIR
bin_PROGRAMS = average_paths epievo_est_complete epievo_est_params_histories \
	epievo_initialization epievo_sim epievo_sim_pairwise \
	extract_state_sequences global_jumps_to_paths scale_tree

average_paths_SOURCES = src/prog/average_paths.cpp
epievo_est_complete_SOURCES = src/prog/epievo_est_complete.cpp
epievo_est_params_histories_SOURCES = src/prog/epievo_est_params_histories.cpp
epievo_initialization_SOURCES = src/prog/epievo_initialization.cpp
epievo_sim_SOURCES = src/prog/epievo_sim.cpp
epievo_sim_pairwise_SOURCES = src/prog/epievo_sim_pairwise.cpp
extract_state_sequences_SOURCES = src/prog/extract_state_sequences.cpp
global_jumps_to_paths_SOURCES = src/prog/global_jumps_to_paths.cpp
scale_tree_SOURCES = src/prog/scale_tree.cpp
