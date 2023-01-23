#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot a histogram of the specified fragment file.

INPUT:  Maximal value for plotting
        Restriction site name
        Species name or description of random sequence, resp.
        Image format ('png' or 'pdf')
        Output file path
        Fragment file
OUTPUT: Image file

@author: Lucia Vedder
@date: 2020, January 21 (Last update: 2021, November 24)
"""

from Bio import SeqIO
from optparse import OptionParser
import matplotlib.pyplot as plt


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
    
    return lengths


def plot_histogram(lengths, name, max_val, rest_name, out_path, f):
    l = list(lengths[name].values())
    l.sort()
    
    v = 0
    for x in l:
        if x < max_val:
            continue
        else:
            v = x
            break
    
    if v == 0:
        v = l[-1]
        
    i = l.index(v)
    
    val_txt = str(max_val/1000) + 'kb'
    bin_size = 100
    if max_val == 500:
        val_txt = '500bp'
        bin_size = 500
    
    title = name + ', cut with ' + rest_name
    out = name.replace('. ', '_') + '_' + rest_name +'_' + val_txt + '.' + f
    plt.hist(l[:i], bin_size, color='green')
    plt.xlabel('Fragment length')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path + out)
    plt.close()



def main():
    # Handle commandline arguments
    parser = OptionParser(usage="%prog [options] <FRAGMENT-FILE>")
    parser.add_option('-m', '--max_value', action='store', dest='max_val', 
                      default=1000, metavar='INT',
                      help="Maximal value for plotting (bp) [default: %default]")
    parser.add_option('-N', '--restriction_name', action='store', dest='rest_name',
                      metavar='STR', default='XX',
                      help="Name of the restriction site enzyme")
    parser.add_option('-f', '--format', action='store', dest='f', 
                      default='png', metavar='FORMAT',
                      help="Image format (one of [\'png\'; \'pdf\']); [default: %default]")
    parser.add_option('-n', '--species_name', action='store', dest='name',
                      metavar='NAME', help="Species name or description of random sequence")
    parser.add_option('-o', '--out_path', action='store', dest='out_path',
                      metavar='PATH', help="Output file path")
    options, args = parser.parse_args()
    
    # Read fragment file and compute lengths
    if options.name:
        name = options.name
    else:
        name = get_name(args[0])
    
    lengths = read_file({name: args[0]})
    
    # Plot histogram
    if options.out_path.endswith('/'):
        out_path = options.out_path
    else:
        out_path = options.out_path + '/'
        
    plot_histogram(lengths, name, int(options.max_val), options.rest_name,
                   out_path, options.f)



if __name__ == '__main__':
    main()
