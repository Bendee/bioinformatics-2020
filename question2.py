#!/usr/bin/python
import time
import sys
import numpy as np

# YOUR FUNCTIONS GO HERE -------------------------------------


SCORE_MAPPING = {
    "A": 3,
    "C": 2,
    "G": 1,
    "T": 2,
    "MISS": -3,
    "-": -4
}


def alignSequences(sequenceOne, sequenceTwo):
    """Locally aligns two given DNA sequences.

    Based on the Needleman-Wunsch algorithm.
    """
    sequences = (sequenceOne, sequenceTwo)
    lengths = [len(sequence) for sequence in sequences]

    # Create and fill scoring and backtracking matrices
    scoreMatrix, backtrackMatrix = fillMatrix(sequences, lengths, *createMatrices(*lengths))

    # Find the starting point for the backracking.
    score, startingPoint = findBestScore(scoreMatrix)
    alignment = backtrackAlignment(sequences, startingPoint, backtrackMatrix)

    return score, alignment


# region Create Matrices
def createMatrices(lengthOne, lengthTwo):
    """ Creates two zero matrices with given dimensions."""
    # Reasoning for methods used:
    #   min(timeit.Timer('[[0 for j in range(10000)] for i in range(10000)]').repeat(5,10))
    #   > 30.831035103999966 seconds
    # vs,
    #   min(timeit.Timer('np.zeros([10000,10000], dtype=int).tolist()').repeat(5,10))
    #   > 14.318026017000193 seconds
    #
    # for the same result.

    # One caveat, lists use more memory than numpy array so I have run into
    # issues running this when I've also had Firefox open.
    scoreMatrix = np.zeros([lengthOne + 1, lengthTwo + 1], dtype=int).tolist()
    backtrackMatrix = np.zeros([lengthOne + 1, lengthTwo + 1], dtype=int).tolist()

    return scoreMatrix, backtrackMatrix


# region Fill Matrices
def fillMatrix(sequences, lengths, scoreMatrix, backtrackMatrix):
    """Calculates scores and fills backtrack matrix."""
    for i in range(1, lengths[0] + 1):
        for j in range(1, lengths[1] + 1):
            # Calculate the potentional scores.
            scores = (
                0,
                scoreMatrix[i - 1][j - 1] + calculateScore(sequences[0][i - 1], sequences[1][j - 1]),
                scoreMatrix[i - 1][j] + SCORE_MAPPING["-"],
                scoreMatrix[i][j - 1] + SCORE_MAPPING["-"]
            )

            # Decide which score to use
            score, direction = decideDirection(scores)
            # The default score and direction in the matrix are 0 
            # and END (0) so we only need to insert when that is
            # not the returned score
            if score != 0:
                scoreMatrix[i][j] = score
                backtrackMatrix[i][j] = direction

    return scoreMatrix, backtrackMatrix


def calculateScore(baseOne, baseTwo):
    """Returns the score for two given genetic bases"""
    if baseOne == baseTwo:
        score = SCORE_MAPPING[baseOne]
    elif baseOne == "-" or baseTwo == "-":
        score = SCORE_MAPPING["-"]
    else:
        score = SCORE_MAPPING["MISS"]

    return score


def decideDirection(scores):
    """Returns the optimal direction for the backtrack matrix."""
    score = max(scores)

    # As index takes the 1st value, it will always choose end or diagonal if it is an option
    # I am using integers to represent the direction in the hope that it will increase speed for initialising matrices
    direction = scores.index(score)

    return score, direction

# endregion Fill Matrices
# endregion Create Matrices


# region Find Best Alignments
def findBestScore(scoreMatrix):
    """Returns the highest score in the scoreMatrix and its position."""
    # Using numpy seemed easiest here, so why not
    scoreMatrix = np.array(scoreMatrix)
    startingPoint = np.unravel_index(
        np.argmax(scoreMatrix, axis=None),
        scoreMatrix.shape
    )
    score = scoreMatrix[startingPoint[0]][startingPoint[1]]

    return score, startingPoint


def backtrackAlignment(sequences, backtrackPoint, backtrackMatrix):
    """Generates an alignment by following the backtrackMatrix from the given point."""
    alignment = ["", ""]

    nextDirection = backtrackMatrix[backtrackPoint[0]][backtrackPoint[1]]
    while nextDirection != 0:
        # 1 is diagonal
        if nextDirection == 1:
            alignment = [
                sequence[backtrackPoint[index] - 1] + alignment[index]
                for index, sequence in enumerate(sequences)
            ]
            backtrackPoint = (backtrackPoint[0] - 1, backtrackPoint[1] - 1)

        # 2 is up
        elif nextDirection == 2:
            alignment = (
                sequences[0][backtrackPoint[0] - 1] + alignment[0],
                "-" + alignment[1]
            )
            backtrackPoint = (backtrackPoint[0] - 1, backtrackPoint[1])

        # 3 is left
        elif nextDirection == 3:
            alignment = (
                "-" + alignment[0],
                sequences[1][backtrackPoint[1] - 1] + alignment[1],
            )
            backtrackPoint = (backtrackPoint[0], backtrackPoint[1] - 1)


        nextDirection = backtrackMatrix[backtrackPoint[0]][backtrackPoint[1]]

    return alignment

# endregion Find Best Alignments


# ------------------------------------------------------------



# DO NOT EDIT ------------------------------------------------
# Given an alignment, which is two strings, display it

def displayAlignment(alignment):
    string1 = alignment[0]
    string2 = alignment[1]
    string3 = ''
    for i in range(min(len(string1),len(string2))):
        if string1[i]==string2[i]:
            string3=string3+"|"
        else:
            string3=string3+" "
    print('Alignment ')
    print('String1: '+string1)
    print('         '+string3)
    print('String2: '+string2+'\n\n')

# ------------------------------------------------------------


# DO NOT EDIT ------------------------------------------------
# This opens the files, loads the sequences and starts the timer
file1 = open(sys.argv[1], 'r')
seq1=file1.read()
file1.close()
file2 = open(sys.argv[2], 'r')
seq2=file2.read()
file2.close()
start = time.time()

#-------------------------------------------------------------


# YOUR CODE GOES HERE ----------------------------------------
# The sequences are contained in the variables seq1 and seq2 from the code above.
# To work with the printing functions below the best local alignment should be called best_alignment and its score should be called best_score. 


best_score, best_alignment = alignSequences(seq1, seq2)


#-------------------------------------------------------------


# DO NOT EDIT ------------------------------------------------
# This calculates the time taken and will print out useful information 
stop = time.time()
time_taken=stop-start

# Print out the best
print('Time taken: '+str(time_taken))
print('Best (score '+str(best_score)+'):')
displayAlignment(best_alignment)

#-------------------------------------------------------------

