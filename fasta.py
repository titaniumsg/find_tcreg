#!/usr/bin/python

#       Sgourakis Lab
#   Author: Santrupti Nerli
#   Modified: Sagar Gupta
#   Date: April 28, 2022
#   Email: snerli@ucsc.edu
#

# import other required libraries
from collections import defaultdict

'''

FASTA class contains all the necessary functionalities required to read
and write to fasta sequence files.

'''

class FASTA:

    # class members
    filename = "" # input file name
    sequences = defaultdict(dict) # sequences read from the fasta file name

    # constructor
    def __init__(self, filename):
        self.filename = filename
        self.sequences = defaultdict(dict)

    # method to get the sequence headers
    def get_headers(self):
        if not self.sequences:
            self.sequences['none'] = ''
        return self.sequences.keys()
        
    # method to get sequences read from fasta file
    def get_sequences(self):
        return self.sequences

    # method to read the fasta file
    # populates headers and sequences
    def read(self):
        if self.filename != None:
            fasta_file_handler = open(self.filename, "r")
            try:
                for line in fasta_file_handler:
                    line = line.rstrip()
                    if ">" in line:
                        header = "".join(list(line)[1:])
                        self.sequences[header] = ''
                    else:
                        self.sequences[header] += line
                fasta_file_handler.close()
            except UnicodeDecodeError:
                print ("UnicodeDecodeError")
                fasta_file_handler.close()
