'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for consensus sequencing

Goes through all reads and creates a counter for reads in a certain offset with most frequent position
Returns new contig list based upon that.

Methods:
1. consensus_sequence()
2. compute_new_contigs

USES the followign dictionary data structure (dictionary within dictionary)
Could imagine it's like a tree
{ CONTIG : {OFFSET : Counter()} }
'''

# IMPROTS
from collections import Counter


'''
Main function to compute new consensus sequence.
Creates data structure object and iterates through each read to 
get a consensus sequence information.

@param reads_dict is a current dictionary with {read : [s, o]}

@return new_contig_sequence
'''
def consensus_sequence(reads_dict):
    # Object
    frequency_info = {}

    # For each read
    for r in reads_dict:
        s = reads_dict[r][0]
        o = reads_dict[r][1]
        # Check if contig number is in structure
        if contig not in frequency_info:
            frequency_info[s] = {}
        for nt in r:
            # Check if offset in the contig number dictionary
            if o not in frequency_info[s]:
                frequency_info[s][o] = Counter()
            # Add counter to that nucleotide in the contig at an offset
            frequency_info[s][o][nt] += 1
            o += 1

    # Send dictionary to compute the new contigs
    return compute_new_contigs(frequency_info)

'''
Compute new contigs.

@param the frequency_info structure with structure { CONTIG : {OFFSET : Counter()} }

@return the new contig sequence
'''
def compute_new_contigs(frequency_info):

    # New list
    new_contigs = []

    # Get sorted list of contigs
    contigs = frequency_info.keys().sort()

    # For each contig
    for s in contigs:
        # Start new contig
        curr = ''
        # Get offsets list
        offsets = frequency_info[s].keys().sort()
        for o in offsets:
            # Get most common nucleotide at certain (s,o)
            nt = frequency_info[s][o].most_common(0)
            curr += nt[0][0]

        new_contigs.append(curr)

    return new_contigs




