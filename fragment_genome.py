#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fragment genomes with one or more specified restriction sites.
If more than one restriction site is given, all combinations are computed.
The output is written to the 'Restriction_fragments' sub-directory.

INPUT:  File path (Genome file has to be in a subfolder called 'Genome')
        List of restriction sites (at least 1)
OUTPUT: FastA file of restriction fragments

@author: Lucia Vedder
@date: May 8, 2019 (Last update: 2021, November 23)
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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


def write_seq(seq, site1, site2, c_name, pos, path):
    if site1 == site2:
        site_path = path + site1 + '/'
        descr = "Cut: " + site1
    else:
        site_path = path + site1 + '-' + site2 +'/'
        descr = "Cut: " + site1 + ", " + site2
    if not os.path.isdir(site_path):
        os.mkdir(site_path)
    
    record = SeqRecord(seq, id=c_name + "[" + str(pos[0]) + ":" + str(pos[1]) + "]",
                        description=descr)
    
    with open(site_path+'all.fa', 'a') as f:
        SeqIO.write(record, f, "fasta")


def fragment(genome, site1, site2, out_path):
    print("Fragment genome:", site1, '-', site2)
    for c_name in genome:
        contig = genome[c_name]
        pos_site1 = []
        pos_site2 = []
        i = 0
        while i < len(contig) and i != -1:
            p_s1 = contig.find(site1, i+1)
            pos_site1.append(p_s1)
            i = p_s1
        j = 0
        while j < len(contig) and j != -1:
            p_s2 = contig.find(site2, j+1)
            pos_site2.append(p_s2)
            j = p_s2
        pos_site1.sort()
        pos_site2.sort()

        start = 0
        if len(pos_site1) < 2 and len(pos_site2) < 2: #no positions for site1 and site2
            end = len(genome[c_name])
            write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
        
        elif len(pos_site1) < 2: #no positions for site1
            end = pos_site2[1]+len(site2)
            write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
            
            j = 1
            pos_site2.append(len(genome[c_name]))
            while j < len(pos_site2)-1:
                start = pos_site2[j]
                end = pos_site2[j+1]+len(site2)
                if end > len(genome[c_name]):
                    end = len(genome[c_name])
                write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
                j = j+1
        
        elif len(pos_site2) < 2: #no positions for site 2
            end = pos_site1[1]+len(site1)
            write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
            
            i = 1
            pos_site1.append(len(genome[c_name]))
            while i < len(pos_site1)-1:
                start = pos_site1[i]
                end = pos_site1[i+1]+len(site1)
                if end > len(genome[c_name]):
                    end = len(genome[c_name])
                write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
                i = i+1
        
        elif site1 == site2: #only one enzyme
            end = pos_site1[1]+len(site1)
            write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
            
            i = 1
            pos_site1.append(len(genome[c_name]))
            while i < len(pos_site1)-1:
                start = pos_site1[i]
                end = pos_site1[i+1]+len(site1)
                if end > len(genome[c_name]):
                    end = len(genome[c_name])
                write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
                i = i+1
        
        else: #2 different enzymes,             
            end = min([pos_site1[1]+len(site1), pos_site2[1]+len(site2)])
            write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)
            
            pos_site2.append(len(genome[c_name]))
            pairs = []
            for i in range(1, len(pos_site1)):
                pairs.append((pos_site1[i], 'site1'))
            for j in range(1, len(pos_site2)):
                pairs.append((pos_site2[j], 'site2'))
            pairs.sort()
            for k in range(len(pairs)-1):
                start = pairs[k][0]
                if pairs[k+1][1] == 'site1':
                    end = pairs[k+1][0]+len(site1)
                else:
                    end = pairs[k+1][0]+len(site2)
                if end > len(genome[c_name]):
                    end = len(genome[c_name])
                write_seq(genome[c_name][start:end], site1, site2, c_name, (start, end), out_path)



def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog <PATH> <LIST,OF,RESTRICTION-SITES>")
    options, args = parser.parse_args()
    
    # Create a list of all genomes in the directory
    paths = get_path_list(args[0])
    print("Number of genomes:", len(paths))
    
    # Run fragmentation for all genomes
    for genome_path in paths:
        genome_path = genome_path.replace('\\', '/')
        genome_file = glob(genome_path+'*')[0]
        handle = SeqIO.parse(genome_file, 'fasta')
        genome = {}
        print("Read genome:", get_name(genome_path))
        for record in handle:
            genome[record.id] = record.seq.upper()
        
        sites = args[1].split(',')
        out_path = genome_path.replace('Genome', 'Restriction_fragments')
        if not os.path.isdir(out_path):
            os.mkdir(out_path)
        
        if len(sites) > 1: #combine all restriction sites
            for site1 in sites:
                for site2 in sites:
                    fragment(genome, site1, site2, out_path)
        else: #just one restriction site
            fragment(genome, sites[0], sites[0], out_path)
    


if __name__ == '__main__':
    main()
