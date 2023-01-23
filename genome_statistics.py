#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute the some basic statistics of the genomes used for the 
in silico restriction fragment analysis.

INPUT:  File path (Genome file has to be in a subfolder called 'Genome')
OUTPUT: -->stdout

@author: Lucia Vedder
@date: 2021, May 26 (Last update: 2021, November 23)
"""

from Bio import SeqIO
from collections import Counter
from optparse import OptionParser
from glob import glob


def get_path_list(path):
    if path.endswith('/'):
        path_list = glob(path+'**/Genome/', recursive=True)
    else:
        path_list = glob(path+'/**/Genome/', recursive=True)
    
    return path_list


def get_name(s):
    parts = s.split('/')
    i = parts.index('Genome')
    return parts[i-1]


def genome_statistics(genome_file):
    handle = SeqIO.parse(genome_file, 'fasta')
    genome = {}
    for record in handle:
        genome[record.id] = record.seq.upper()
    
    seq = ""
    lengths = []
    num = 0
    for contig in genome:
        seq = seq + str(genome[contig])
        num += 1
        lengths.append(len(genome[contig]))
    
    base_counts = Counter(seq)
    gc = (base_counts['C'] + base_counts['G']) / len(seq) * 100
    
    lengths.sort()
    n = 0
    n50 = 0
    for l in lengths:
        n += l
        if n >= len(seq)/2:
            n50 = l
            break
    
    print("Genome size:", len(seq))
    print("GC-content:", gc, "%")
    print("Number of contigs:", num)
    print("N50:", n50)
    
    
def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog <PATH>")
    options, args = parser.parse_args()
    
    # Create a list of all genomes in the directory
    paths = get_path_list(args[0])
    print("Number of genomes:", len(paths))
    
    # Run statistics for all genomes
    for genome_path in paths:
        genome_path = genome_path.replace('\\', '/')
        genome_file = glob(genome_path+'*')[0]
        
        print("\n=====", get_name(genome_path), "=====")
        genome_statistics(genome_file)



if __name__ == '__main__':
    main()
