#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Search for peaks in the restriction fragment distributions.
A peak is defined as a count from single fragment length that deviates more 
than +-2*stdev from the mean of all counts up to 500bp.

INPUT:  Species name or description of random sequence, resp.
        Fragment file
OUTPUT: --> stdout (peak info)
        FastA files

@author: Lucia Vedder
@date: 2021, August 03 (Last update: 2021, November 24)
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import Counter
from optparse import OptionParser
import numpy as np


def get_name(s):
    parts = s.split('/')
    i = parts.index('Restriction_fragments')
    return parts[i-1]


def read_file(files):
    fragments = {}
    for name in files:
        handle = SeqIO.parse(files[name], 'fasta')
        f = {}
        for record in handle:
            f[record.id] = record.seq.upper()
        fragments[name] = f
    
    lengths = {}
    for name in fragments:
        l = {}
        for n in fragments[name]:
            l[n] = len(fragments[name][n])
        lengths[name] = l
    
    return lengths, fragments


def find_peaks(lengths, fragments, files):
    for name in lengths:
        counts = Counter(lengths[name].values()).items()
        selected = []
        for pair in counts:
            if pair[0] <= 500:
                selected.append(pair)
        data = dict(selected)
        
        data_mean = np.mean(list(data.values()))
        data_stdev = np.std(list(data.values()))
        sig_keys = []
        for key in data:
            if data[key] > data_mean + 2 * data_stdev:
                sig_keys.append(key)
        
        # Print peak info
        print("\n=====", name, "=====")
        print("Peaks: (fragment length: count)")
        for key in sig_keys:
            print('(', key, ': ', data[key], ')', sep='')
        
        # Write peak sequences
        path = files[name].replace('all.fa', '')
        for key in sig_keys:
            file_name = 'peaks_length-' + str(key) + '.fa'
            
            for n in lengths[name]:
                if lengths[name][n] == key:
                    record = SeqRecord(fragments[name][n], id=n)
                    with open(path + file_name, 'a') as f:
                        SeqIO.write(record, f, 'fasta')



def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog [options] <FRAGMENT-FILE>")
    parser.add_option('-n', '--species_name', action='store', dest='name',
                      metavar='NAME', help="Species name or description of random sequence")
    options, args = parser.parse_args()
    
    # Read fragment file and compute lengths
    if options.name:
        name = options.name
    else:
        name = get_name(args[0])
    
    lengths, fragments = read_file({name: args[0]})
    
    # Count fragment lengths and compute peaks
    # Write info to stdout and save sequences as fasta
    find_peaks(lengths, fragments, {name: args[0]})



if __name__ == '__main__':
    main()
