'''
Dan Adler, Farhan Damani
Computational Genomics Final

Library for contig merging.  Contains 3 methods:

1. merge_check_global() which goes through all n choose 2 contigs and checks whether a merge should occur
2. merge_check_local() checks two contigs to see if a merge should occur
3. merge_contigs() which performs an actual merge of contigs necessary to be merged

NEED: To somehow look at merged "paths" (meaning more than 2 contigs merged at once)...unless we can make
      make this not happen

NEED: To look at if there is a contig that aligns to two other contigs what happens

IDEA for this.  Keep track of contigs to be merged in a list.  If it appears in a list that it will be merged
look at how it would be merged twice.  I think I need like a merge 'prep' function step, where I look at how the merges
would occur.  Think I'm going to write this function with a test case.
'''

# IMPORTS
import itertools

# Global variables

'''
Checks all nchoosek(contig_length,2) contigs using merge_check_local() to see which contigs
should be merged using merge_contigs()

@param contigs_list is the list of all current contigs
@param OVERLAP_LENGTH (defaults to 15)

@return the overlaps or a blank list
'''
def merge_check_global(contig_list, OVERLAP_LENGTH = 15):
    # Merge will be 0, and a 1 is added each time merge_local() comes back as a needed merge
    merge = 0
    # Dictionary to hold permutations
    contig_dict = {}
    # Go through each contig
    for suffix in contig_list:
        # Check if a merge is needed using merge_local()
        contig_dict[suffix] = Counter()
        for prefix in conting_list:
            if prefix != suffix:
                contig_dict[suffix][prefix] = suffixPrefixMatch(suffix, prefix, OVERLAP_LENGTH)
        merge += contig_dict[list(pair)]

    # Check if merge_contigs needs to occur
    if merge > 0:
        return contig_dict
    else:
        return []


'''
Finds the longest suffix that is a prefix of
another string.
@param str1 is the string whose suffix will be checked
@param str2 is the string whose prefix will be checked
@param min_overlap is the minimum overlap length
@param OVERLAP_LENGTH (defaults to 15)

@returns the length of the match
'''
def suffixPrefixMatch(str1, str2, min_overlap):
    ''' Returns length of longest suffix of str1 that is prefix of
        str2, as long as that suffix is at least as long as min_overlap. '''
    if len(str2) < min_overlap: return 0
    str2_prefix = str2[:min_overlap]
    str1_pos = -1
    while True:
        str1_pos = str1.find(str2_prefix, str1_pos + 1)
        if str1_pos == -1: return 0
        str1_suffix = str1[str1_pos:]
        if str2.startswith(str1_suffix): return len(str1_suffix)

'''
Suffix filter (i.e. choose greatest prefix match for a suffix)

@param contigs_dict is a dictionary for each contig as a suffix that has a counter object
@return each contig in a dictionary with their respective right best buddy
'''
def suffix_filter(contig_dictionary):
    bbr = {}

    for suffix in contig_dictionary:
        if len(contig_dictionary[suffix]) > 0:
            c = contig_dictionary[suffix].most_common(1)
            bbr[c] = [c[[0][0], c[0][1]]

    return bbr

'''
Compute ordering of new  contigs.

@param bbr is the best budy dicitonary
@return a list of triplets that show [(contig, next, overlap)...]
'''
def compute_order(bbr):
    # Somehow take walk through things and find orders
    pass

'''
Perform a merge contigs.
Go through the entered dictionary with a 1 if a contigs needs to merged and 2 otherwise

@param contig_list is the list of contigs
@param contig_dictionary is a dicitonary with key = [contig1,contig2] and value
@param reads_dict is a dictionary of each read as: {read: [[s_0,o_0]...}

@return [new_contig_list, new_reads_dict] is the new list of contgs
'''
def merge_contigs(contig_list, contig_dictionary, reads_dict, OVERLAP_LENGTH = 15):
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
            # Change reads dict
            new_reads_dict = change_reads_on_merge(reads_dict, pair[0], pair[1], len(new_contig_list) - 1, len(contig1))

    return [new_contig_list, new_reads_dict]

'''
Change position of reads during a merge

@param reads_dict is the dictionary of reads (assuming {read: [s, o]})
@param contig1 is the offset for contig1 merged
@param contig2 is the offset for contig2 merged
@param new_s is the new offset to use for the contig (updated s)
@param contig1_length is the length of contig1 (to compute new offsets)

@return reads_dict once updated
'''
def change_reads_on_merge(reads_dict, contig1, contig2, new_s, contig1_length):

    for read in reads_dict:
        for i in range(len(reads_dict[read])):
        # If in first contig just need to change contig_pos (s)
        if reads_dict[read][i][0] == contig1:
            reads_dict[read][i][0] = new_s
        # If in second contig need to change s and o
        elif reads_dict[read][i][0] == contig2:
            reads_dict[read][i][0] = new_s
            reads_dict[read][i][1] = contig1_length + reads_dict[read][1] - OVERLAP_LENGTH - 1

    return reads_dict



