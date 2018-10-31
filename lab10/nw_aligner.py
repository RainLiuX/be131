#!/usr/bin/python

"""

Needleman-Wunsch Aligner
Bioengineering 131/231, Fall 2018

Command-line script to read in the contents of a multi-FASTA containing two
sequences and a score matrix, perform a global alignment, and print the
resulting alignment matrix and optimal alignment to STDOUT.

"""

import os
import sys
from Bio import SeqIO

class NWAligner:
    def __init__(self, score_matrix_fname):
        self.score_matrix, self.gap_penalty = self.load_score_matrix(score_matrix_fname)

    @staticmethod
    def load_score_matrix(fname):
        """
        Input: (String) A path to a scoring matrix file.
        Output: (Dictionary) A nested dictionary where keys are strings
                and elements are scores as integers.
    
        Example:
    
        >>> matrix, gap_penalty = NWAligner.load_score_matrix('/home/bioe131/BLOSUM62')
        >>> matrix['A']['A']
        4
        >>> matrix['W']['W']
        11
        >>> gap_penalty
        -4

        """

        score_matrix = {}
        gap_penalty = None
        allowed_letters = []

        with open(fname) as fp:
            for line_num, line in enumerate(fp):
                # ignore comments in matrix file
                if line.startswith("#"):
                    continue
                # read allowed letters from first line
                if line_num == 0 :
                    allowed_letters = line.rstrip().split(' ')
                    # initialize sub dictionaries
                    for letter in allowed_letters :
                        score_matrix[letter] = {}
                # read gap penalty from last line
                elif line_num == len(allowed_letters)+1 :
                    gap_penalty = int(line.rstrip())
                else :
                    # read score matrix
                    line_split = line.rstrip().split(' ')
                    for i in range(line_num) :
                        # force symmetric
                        score_matrix[allowed_letters[line_num-1]][allowed_letters[i]] = int(line_split[i])
                        score_matrix[allowed_letters[i]][allowed_letters[line_num-1]] = int(line_split[i])
        return score_matrix, gap_penalty

    @staticmethod
    def load_FASTA(fname):
        """
        Input: (String) A path to a FASTA file containing exactly two sequences.
        Output: (List) A list containing two strings: one for each sequence.

        Example:

        >>> seqs = NWAligner.load_FASTA('example.fa')
        >>> seqs[0]
        'YAADSKATPGNPAFHQDEIFLARIAFIYQMWDGGQLKLIDYAPHHVMCEE'
        >>> seqs[1]
        'WVGQPNMKVQHWSNMKACCVKFITWTFIAPEKHACKWTETAYQADCDIIW'
        >>> len(seqs)
        2

        """

        seqs = []

        for rec in SeqIO.parse(fname, 'fasta') :
            seqs.append(rec.seq)
        if len(seqs) > 2 :
            sys.exit(1)

        return seqs

    def align(self, seq_x, seq_y, print_matrix = False):
        """
        Input: (Strings) Two sequences to be aligned (seq_x and seq_y).
               (Boolean) If print_matrix is True, print the dynamic programming
                         matrix before traceback.
        Output: (Tuple of strings) Two sequences, aligned.

        Example:

        >>> aligner = NWAligner('BLOSUM62')
        >>> seqs = aligner.load_FASTA('example.fa')
        >>> aligner.align(seqs[0], seqs[1])
        ('YAAD-SKATPGNPAF---HQDEIF--L-AR--IA-FIYQM-WDGGQLK-LIDYAPH-HVM-C---E-------E---',
         'W---VGQ--P-N--MKVQH----WSNMKA-CCV-KFI---TW------TFI--APEKH--ACKWTETAYQADCDIIW')

        """

        ###
        ### INITIALIZATION
        ###

        # create two empty matrices with sizes based on the input sequences.
        # one contains the dynamic programming matrix, the other contains
        # pointers we'll use during traceback
        matrix = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]
        pointers = [[0] * (len(seq_y) + 1) for _ in range(len(seq_x) + 1)]

        # fill in the first column(seq_x) and first row(seq_y)
        for i in range(len(seq_y)+1) :
            matrix[0][i] = self.gap_penalty * i
        for i in range(len(seq_x)+1) :
            matrix[i][0] = self.gap_penalty * i

        ###
        ### RECURSION
        ###

        # fill the dynamic programming and pointer matrices
        for x in range(1, len(seq_x) + 1):
            for y in range(1, len(seq_y) + 1):
                match_score = self.score_matrix[seq_x[x - 1]][seq_y[y - 1]]
                score = max(\
                        matrix[x - 1][y - 1] + match_score,\
                        matrix[x - 1][y] + self.gap_penalty,\
                        matrix[x][y - 1] + self.gap_penalty)
                matrix[x][y] = score

                if score == matrix[x - 1][y] + self.gap_penalty :
                    pointers[x][y] = 1 # 1 goes down
                elif score == matrix[x - 1][y - 1] + match_score :
                    pointers[x][y] = 2 # 2 goes diagonal
                elif score == matrix[x][y - 1] + self.gap_penalty :
                    pointers[x][y] = 3 # 3 goes right

        # print the dynamic programming matrix
        if print_matrix:
            for x in range(len(seq_x) + 1):
                # following line is edited for py3
                print(" ".join(map(lambda i: str(int(i)), matrix[x])))

        ###
        ### TRACEBACK
        ###

        # starting from the bottom right corner, follow the pointers back
        x, y = len(seq_x), len(seq_y)

        # fill these lists with the aligned sequences
        align_x = []
        align_y = []

        while x > 0 or y > 0:
            move = pointers[x][y]
            if move == 0 :
                # move to the edge
                break
            elif move == 1:
                # down, x is aligned to gap
                align_x.append(seq_x[x - 1])
                align_y.append('-')
                x -= 1
            elif move == 2:
                # diagonal, x is aligned to y
                align_x.append(seq_x[x - 1])
                align_y.append(seq_y[y - 1])
                x -= 1
                y -= 1
            elif move == 3:
                # right, y is aligned to gap
                align_x.append('-')
                align_y.append(seq_y[y - 1])
                y -= 1

        # flip the alignments, as they're reversed
        return ("".join(align_x[::-1]), "".join(align_y[::-1]))

###                                      ###
### NO NEED TO EDIT CODE BELOW THIS LINE ###
###                                      ###

if __name__ == '__main__':
    def usage():
        print('usage: %s matrixfilename stringfilename')
        sys.exit(1)

    if len(sys.argv) != 3:
        usage()

    for fname in sys.argv[1:]:
        if not os.path.isfile(fname):
            print('Can not open %s' % (fname,))
            usage()

    aligner = NWAligner(sys.argv[1])
    seqs = aligner.load_FASTA(sys.argv[2])
    result = aligner.align(seqs[0], seqs[1])

    print('>seq1\n%s\n>seq2\n%s' % (result[0], result[1]))
