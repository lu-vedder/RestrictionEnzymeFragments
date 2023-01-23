#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute statistical values for the restriction fragment distributions.

INPUT:  Restriction fragment file path
OUTPUT: --> stdout

@author: Lucia Vedder
@date: 2021, May 18 (Last update: 2021, November 23)
"""

from Bio import SeqIO
from statistics import mean, stdev, mode
from optparse import OptionParser
from glob import glob


def get_path_list(path, label):
    if path.endswith('/'):
        path_list = glob(path+'**/'+label+'/', recursive=True)
    else:
        path_list = glob(path+'/**/'+label+'/', recursive=True)
    
    return path_list


def get_name(s, label):
    parts = s.split('/')
    i = parts.index(label)
    return parts[i-1]


def fragment_statistics(fragment_file):
    handle = SeqIO.parse(fragment_file, 'fasta')
    fragments = {}
    for record in handle:
        fragments[record.id] = record.seq
    
    lengths = {}
    for name in fragments:
        lengths[name] = len(fragments[name])
        
    print("Number of fragments:", len(fragments))
    print("Min. fragment length:", min(lengths.values()))
    print("Max. fragment length:", max(lengths.values()))
    print("Average fragment length:", round(mean(lengths.values()), 3))
    print("Standard deviation:", round(stdev(lengths.values()), 3))
    print("Most frequent fragment length:", mode(lengths.values()))
    
    return list(lengths.values())



def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog <PATH>")
    options, args = parser.parse_args()
    
    # Create a list of all fragmented genomes in the directory
    paths = get_path_list(args[0], 'Restriction_fragments')
    print("Number of fragmented genome files:", len(paths))
    
    # Run statistics for all fragment files
    for p in paths:
        p = p.replace('\\', '/')
        fragment_path_list = get_path_list(p, '*') #get all restriction sites
        
        for fragment_path in fragment_path_list:
            fragment_file = fragment_path.replace('\\', '/') + 'all.fa'
                    
            print("\n=====", get_name(p, 'Restriction_fragments'), "-", get_name(fragment_file, 'all.fa'), "=====")
            fragment_statistics(fragment_file)


if __name__ == '__main__':
    main()
