'''
Dan Adler, Farhan Damani
Computational Genomics Final

CreateTestOutput.py

COMMAND:
python CreateInputData.py n k input.fa output.fa n

Takes a genome (all caps no spaces) and
randomly samples n reads of length k from that genome.
'''

# Imports
import sys
import random

'''
Opens file into a string.
@param infile is the file to be opened
@return the output string (string from the infile)
'''
def openFile(infile):
    data = open(infile, 'r').read()
    return data

'''
Samples n random reads of length k from that genome.
@param data is the input genome.
@param n is the number of sample reads.
@param k is the length of reads to sample.
@return reads of that data.
'''
def sampleReads(data, n, k):
    reads = []
    for i in range(n):
        start = random.randint(0,len(data)-100)
        reads.append('>' + str(i+1))
        reads.append(data[start:start+100])
    return '\n'.join(reads)

'''
Output string we want to file
@param outfile name is the name of the output file
@param output_string is the otuput string
'''
def outputFile(outfile, output_string):
    of = open(outfile, 'w')
    of.write(output_string)
    of.close()

'''
Main function
'''
def main():
    # Open file
    data = openFile(sys.argv[3])
    # Get random reads
    reads = sampleReads(data,int(sys.argv[1]),int(sys.argv[2]))
    # Output file
    outputFile(sys.argv[4],reads)

main()