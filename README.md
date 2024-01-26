# RestrictionEnzymeFragments

These python scripts are used for the restriction fragment analysis pubished by Vedder & Schoof (2024) [1].


**genome_statistics.py**  
Compute the some basic statistics of the genomes used for the in silico restriction fragment analysis.

**random_sequence.py**  
Create a random sequence with a specified GC-content (default 0.5).

**fragment_genome.py**  
Fragment genomes with one or more specified restriction sites. If more than one restriction site is given, all combinations are computed. The output is written to the 'Restriction_fragments' sub-directory.

**fragment_statistics.py**  
Compute statistical values for the restriction fragment distributions.

**plot_fragments.py**  
Plot a histogram of the specified fragment file.

**find_peaks.py**  
Search for peaks in the restriction fragment distributions. A peak is defined as a count of single fragment length that deviates more than +-2*stdev from the mean of all counts up to 500bp.

# Citation
Please cite via Zenodo: [![DOI](https://zenodo.org/badge/563350631.svg)](https://zenodo.org/doi/10.5281/zenodo.10572873)

# License
Copyright (c) 2024 Lucia Vedder <br>
For details see the [LICENSE](LICENSE) file.

---
[1] Vedder, L. and Schoof, H. (2024): In silico restriction site analysis of whole genome sequences shows patterns caused by selection and sequence duplications. Heliyon. In review.
