#!/usr/bin/python

import os
import sys

header = "#E      C       D"

rootseq = sys.argv[1]
output = sys.argv[2]

OUT = open(output, "w")

print >>OUT, header
for (pos, state) in enumerate(rootseq):
    print >>OUT, "\t".join([str(pos), state, '0', '0'])

OUT.close()
