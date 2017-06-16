# File: NeedlemanWunsch.py
# Author: Jason Lu (jasonlu6@bu.edu)
# Collaborations:
# Rohim Banerji (rohinb96@bu.edu) (Class of 2019 Fellow Undergraduate CE student)
# Professor Douglas Densmore (dougd@bu.edu)
# Marliene Pavan (mapavan@bu.edu) (Supervisor / Senior Researcher at BU)
# Densmore CIDARLAB UROP Project #2 (WETLAB division)
# Date: 6/13/2017

# This program is an application of Needleman Wunsch Algorithm (with
# actual genetic sequences from Mary's DAMP Lab data to determine which
# one of the sequences are the closest from the chosen subset of 10 sequences)

# Recursive implementation for the helper functions, iterative using Dynamic Programming
# for the main NW function

# Bug to fix: List Index out of range error

############################################################################
# imports:

# make sure numpy is within the function!
import numpy as np
from string import *

# helper function to determine case match or mismatch
# parameters:
# n1: sequence 1 of DNA (by characters)
# n2: sequence 2 of DNA (by characters)
# point: the point of intersection between the n1,n2 (as a list)

def MatchOrMismatch(n1, n2, point):
    point = ['MATCH','MISMATCH',' ']
    # it's a match!
    if (n1 == n2):
        print "It's a match!"
        return point[0]
    # it's a mismatch!
    elif (n1 != n2):
        print "It's a mismatch!"
        return point[1]
    # else, undetermined
    else:
        print "Sequence part not determined"
        return point[2]

# helper function to give optimal elements of the algnemnt (F-matrix) grid, and
# returns pointers (V = '|', H = '-', D = '/')  for the pointers matrix

def PointerDiagram(diagonal,vertical,horizontal):

    # maximum absolute distance (edit-distance) between each pointer
    pointer = max(diagonal,vertical,horizontal)

    # diagonal is the max
    if (diagonal  == pointer):
        return '/'

    # horizontal is the max
    elif (horizontal == pointer):
        return '-'

    # vertical is the max
    elif (vertical == pointer):
        return '|'

    # otherwise, return nothing
    else:
        return ' '

# main function that uses Dynamic Programming for the Needleman-Wunsch algorithm

# Penalties:
# match += 1 (each match will reduce the inaccuracy by 1 full point
# mismatch += -1 (each mismatch will reduce the inaccuracy by 1 full point)
# gap/missing char += -2 (worst case, each gap will reduce the inaccuracy by 2 full points)
def NeedlemanWunsch(seq1,seq2,match=1,mismatch=-1,gap=-2):

    # dictionary to store the penalites:
    penalty = {'MATCH': match, 'MISMATCH': mismatch, 'GAP': gap}
    # dimension of matrix
    # columns: avoid one-off error
    col = len(seq1) + 1
    # rows: avoid one-off error
    row = len(seq2) + 1
    # init alignment matrix with 0's
    align_matrix = np.zeros((row,col),dtype = int)
    # init pointer matrix with '0' (as a string)
    pointer_matrix = np.zeros((row,col),dtype = str)

    # for loop to go through first rows element in the matrix, and assign the gap penalty
    for r in range(row):
        align_matrix[r][0] = penalty['GAP'] * r
        pointer_matrix[r][0] = '|'

    # for loop to go through the first columns element in the matrix, assign the gap penalty
    for c in range(col):
        align_matrix[0][c] = penalty['GAP'] * c
        pointer_matrix[0][c] = '-'

    # filling the matrix with the recursive formula approach:

    # Recurrence relation:

    # Parameters:
    # Matrix: the initial and subsequent grid
    #

    # Matrix(row,col) = max {Matrix(row-1,col+1) + S(A[row],B[col]),
    # Matrix(row,col-1) + dist, Matrix(row-1,col) + dist)

    # return the first element of pointer element matrix
    pointer_matrix[0][0] = 0
    print pointer_matrix
    print align_matrix

    # for loop to go through the rows:
    for rows in range(row):
        # for loop to go through columns:
        for cols in range(col):
            # recursively build up / memoize the diagonals for mismatch / match
            print penalty
            # diag = align_matrix[rows-1][cols-1] + MatchOrMismatch(seq1[cols-1],seq2[rows-1],str(penalty['MISMATCH']))
            # print "The value for diagonal match / mismatch: " , diag
            # build up the value for gap (horizontal)
            horiz = align_matrix[rows][cols-1] + penalty['GAP']
            print "The value of gap (horizontal): " , horiz
            # build up the value for gap (vertical)
            vert = align_matrix[rows-1][cols] + penalty['GAP']
            print "The value of gap (vertical): " , vert
            #align_matrix[rows][cols] = max(diag,horiz,vert)
            #pointer_matrix[rows][cols] = PointerDiagram(diag,vert,horiz)

    # print out the alignment matrix:
    print align_matrix
    print pointer_matrix
    # np.matrix(align_matrix)

    # print out of the pointer matrix:
    print np.matrix(pointer_matrix)

    # main function to test

def main():

    # first test
    sequence1 = "TGCTAGCA"
    sequence2 = "TCCCGGATA"
    sequence3 = "ACTACTACTA"
    sequence4 = "ACTGACTGACTG"

    print "Length of sequence 1: ", len(sequence1)
    print "Length of sequence 2: ", len(sequence2)
    print "Length of sequence 3: ", len(sequence3)
    print "Length of sequence 4: ", len(sequence4)

    maxSeq = max(len(sequence1),len(sequence2))

    # something is wrong with this function?!
    print MatchOrMismatch(sequence1, sequence2,[maxSeq])

    # print out the diagram:
    print PointerDiagram(diagonal=10,vertical=10,horizontal=10)

    # second test
    print NeedlemanWunsch(sequence1,sequence2,match=1,mismatch=-1,gap=-2)

main()












