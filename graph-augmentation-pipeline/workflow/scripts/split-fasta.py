#!/usr/bin/env python

import sys
import re

not_first = False

def process_seq(header, seq):
    # reduce repeated Ns                                                                                                                                                                                                                          
    seq = re.sub(r'[^\w\s]|(.)(?=\1)', '', seq)
    # get indices of Ns                                                                                                                                                                                                                           
    match_idx = 1
    posn_idx = 0
    for match in re.finditer(r'N(.*?)', seq):
        if seq[posn_idx:match.start()]:
        	sys.stdout.write("%s_%d\n%s\n" % (header, match_idx, seq[posn_idx:match.start()]))
        posn_idx = match.end()
        match_idx += 1
    if seq[posn_idx:]:
        sys.stdout.write("%s_%d\n%s\n" % (header, match_idx, seq[posn_idx:]))

for line in sys.stdin:
    if line.startswith('>'):
        if not_first:
            process_seq(header, seq)
        header = line.rstrip()
        seq = ''
        not_first = True
    else:
        seq += line.rstrip()

if not_first:
    process_seq(header, seq)
