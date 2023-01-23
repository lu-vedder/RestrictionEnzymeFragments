#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a random sequence with a specified GC-content (default 0.5).

INPUT:  Sequence length
        GC-content [0 - 1]
        Output file name
OUTPUT: FastA file

@author: Lucia Vedder
@date: 2020, January 21 (Last update: 2021, November 24)
"""

import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from optparse import OptionParser


def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog [options]")
    parser.add_option('-l', '--length', action='store', dest='seq_len',
                      help="Length of the random sequence (bp) [default: %default]",
                      default=1000, metavar='INT')
    parser.add_option('-o', '--out', action='store', dest='out', metavar='FILE',
                      help="Output file (fasta)")
    parser.add_option('-c', '--gc-content', action='store', dest='gc_percent', 
                      help="GC-content [0 - 1]; [default: %default]", 
                      default=0.5, metavar='NUM')
    options, args = parser.parse_args()
    
    # Create random sequence
    gc = 50 * float(options.gc_percent)
    at = 50 - gc
    record = SeqRecord(Seq(''.join([random.choices(['A','C','G','T'], weights=[at, gc, gc, at])[0] for x in range(options.seq_len)])), id="random")    
    with open(options.out, 'w') as f:
        SeqIO.write(record, f, "fasta")



if __name__ == '__main__':
    main()
