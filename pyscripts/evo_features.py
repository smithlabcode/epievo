#!/usr/bin/python
#
# Copyright (C) 2019-2020 University of Southern California and Liz Ji
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import argparse

#=================================================================
def read_states(FIN):
    block_states = []
    block_sizes = []

    state = 4 # flag of starting
    prev_state = state
    size = 0

    FIN.readline()    # skip header line

    for line in FIN:
        (pos, root, leaf) = [int(i) for i in line.strip().split()]
        state = 2 * root + leaf

        if not state == prev_state:    # new block
            if not prev_state == 4:    # going to store last block
                block_sizes.append(size)
                block_states.append(prev_state)
            size = 1
        else:
            size += 1

        prev_state = state

    # store last block
    if not prev_state == 4:
        block_sizes.append(size)
        block_states.append(prev_state)

    return (block_states, block_sizes)


def assign_features(block_states):
    def check_triple(l, m, r):
        feature = 'X'
        if m == 0:
            feature = 'N'
        elif m == 3:
            feature = 'F'
        elif m == 1:
            if l < 3 and r < 3:
                feature = 'B'
            elif l == 3 and r == 3:
                feature = 'M'
            else:
                feature = 'E'
        else:    # m = 2
            if l < 3 and r < 3:
                feature = 'D'
            elif l == 3 and r == 3:
                feature = 'S'
            else:
                feature = 'C'
        return feature

    def check_pair(m, lr):
        feature = 'X'
        if m == 0:
            feature = 'N'
        elif m == 3:
            feature = 'F'
        elif m == 1:
            if lr < 3:
                feature = 'B'
            else:
                feature = 'E'
        else:    # m = 2
            if lr < 3:
                feature = 'D'
            else:
                feature = 'C'
        return feature

    def check_single(m):
        feature = 'X'
        if m == 0:
            feature = 'N'
        elif m == 1:
            feature = 'B'
        elif m == 2:
            feature = 'D'
        else:
            feature = 'F'
        return feature

    #######################################
    block_features = []
    if len(block_states) == 1: # only one block
        block_features.append(check_single(block_states[0]))
    else:
        block_features.append(check_pair(block_states[0],block_states[1]))
        if len(block_states) > 2: # at least three blocks
            for i in xrange(1, len(block_states)-1):
                block_features.append(check_triple(block_states[i-1],
                                                   block_states[i],
                                                   block_states[i+1]))
        block_features.append(check_pair(block_states[-1],block_states[-2]))
    return block_features

#=================================================================
#=================================================================
#=================================================================
#=================================================================

def main():

    #=============================================================
    parser = argparse.ArgumentParser(description = 'counting\
        evolution features')
    parser.add_argument('-i', '--input', dest='input', required=True,\
                        help='input file (.outmeth)')

    args = parser.parse_args()

    #=============================================================

    # (1). Store block sizes and states
    # 0: o    1: o    2: x    3: x
    #    o       x       o       x
    FIN = open(args.input)
    (block_states, block_sizes) = read_states(FIN)

    # (2). Assign block feature
    # F: Fill        N: None
    # B: Birth       D: Death
    # E: Expansion   C: Contraction
    # M: Merge       S: Split
    block_features = assign_features(block_states)

    # (3). Count features and output statistics
    features = ['B', 'D', 'E', 'C', 'M', 'S']
    to_report = ["A_00", "A_01", "A_10", "A_11"] + \
                ["N_"+f for f in features] + \
                ["A_"+f for f in features]
    report = dict(zip(to_report, [0.0]*len(to_report)))

    for i in xrange(0, len(block_states)):
        if block_states[i] == 0:
            report["A_00"] += block_sizes[i]
        elif block_states[i] == 1:
            report["A_01"] += block_sizes[i]
        elif block_states[i] == 2:
            report["A_10"] += block_sizes[i]
        else:
            report["A_11"] += block_sizes[i]
        if not block_features[i] in ['N', 'F']:
            report["N_"+block_features[i]] += 1
            report["A_"+block_features[i]] += block_sizes[i]

    print "\t".join(to_report) + "\n" + \
          "\t".join([str(report[t]) for t in to_report])

if __name__ == "__main__":
    main()
