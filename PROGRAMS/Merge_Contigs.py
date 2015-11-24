'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for contig merging.  Contains 3 methods:

1. merge_check_global() which goes through all n choose 2 contigs and checks whether a merge should occur
2. merge_check_local() checks two contigs to see if a merge should occur
3. merge_contigs() which performs an actual merge of contigs necessary to be merged


NEED: To somehow updates where the reads are positioned within the contigs once two contigs are merged

'''

# IMPORTS
import itertools

# Global variables
OVERLAP_LENGTH = 15

'''
Checks all nchoosek(contig_length,2) contigs using merge_check_local() to see which contigs
should be merged using merge_contigs()

@param contigs_list is the list of all current contigs
@return new_contigs_list as the new contigs list.
'''
def merge_check_global(contig_list):
    # Merge will be 0, and a 1 is added each time merge_local() comes back as a needed merge
    merge = 0
    # Dictionary to hold permutations
    contig_dict = {}
    # Generate all nchoose2 pairs of numbers between 1 through length of contigs
    # Go through each pair
    for pair in itertools.permutations(range(len(contig_list)),2):
        # Check if a merge is needed using merge_local()
        contig_dict[list(pair)] = merge_check_local(contig_list[pair[0]], contig_list[pair[1]])
        merge += contig_dict[list(pair)]

    # Check if merge_contigs needs to occur
    if merge > 0:
        return merge_contigs(contig_list)
    else:
        return contig_list

'''
Check if two contigs need to be merged (i.e. overlap > k) where k = 15 as a start.

@param contig1 is the contig whose suffix will be checked
@param contig2 is the contig whose prefix will be checked
@return 0 if overlap is not correct and 1 if overlap is a correct alignment
'''
def merge_check_local(contig1, contig2):
    # Check if the OVERLAP_LENGTH parameters are the same
    if contig1[-OVERLAP_LENGTH:] == contig2[0:OVERLAP_LENGTH]:
        return 1
    else:
        return 0

'''
Perform a merge contigs.
Go through the entered dictionary with a 1 if a contigs needs to merged and 2 otherwise

@param contig_list is the list of contigs
@param contig_dictionary is a dicitonary with key = [contig1,contig2] and value
@return new_contig_list is the new list of contgs
'''
def merge_contigs(contig_list, contig_dictionary):
    # Create new list
    new_contig_list = []
    # Go through each pair
    for pair in contig_dictionary:
        # If a merge needs to take place between two contigs
        if contig_dictionary[pair] == 1:
            # Get first contig
            contig1 = contig_list[pair[0]]
            # Get second contig
            contig2 = contig_list[pair[1]]
            # Merge with overlap as OVERLAP_LENGTH
            new_contig_list.append(contig1 + contig2[OVERLAP_LENGTH:])

    return new_contig_list


