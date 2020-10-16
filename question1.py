#!/usr/bin/python
import time
import sys


# YOUR FUNCTIONS GO HERE -------------------------------------


SCORE_MAPPING = {
    "A": 3,
    "C": 2,
    "G": 1,
    "T": 2,
    "MISS": -3,
    "-": -4
}
num_alignments = 0


def alignSequences(sequenceOne, sequenceTwo):
    """Globally aligns to DNA sequences."""
    global num_alignments

    # Check if the end points have been reached
    if sequenceOne == "":
        length = len(sequenceTwo)
        score = length * calculateScore("-", "-")
        alignment = ("-"*length, sequenceTwo)
        num_alignments += 1
    elif sequenceTwo == "":
        length = len(sequenceOne)
        score = length * calculateScore("-", "-")
        alignment = ("-"*length, sequenceOne)
        num_alignments += 1

    else:
        # Calculate the possible scores and alignments
        alignScoreOne, alignmentOne = alignSequences(sequenceOne[:-1], sequenceTwo[:-1])
        alignScoreTwo, alignmentTwo = alignSequences(sequenceOne[:-1], sequenceTwo)
        alignScoreThree, alignmentThree = alignSequences(sequenceOne, sequenceTwo[:-1])
        potentialAlignments = (alignmentOne, alignmentTwo, alignmentThree)

        scores = (
            alignScoreOne + calculateScore(sequenceOne[-1], sequenceTwo[-1]),
            alignScoreTwo + calculateScore("-","-"),
            alignScoreThree + calculateScore("-","-")
        )
        score = max(scores)
        index = scores.index(score)

        # Form the alignment to return
        alignment = potentialAlignments[index]
        if index == 0:
            alignment = (
                alignment[0] + sequenceOne[-1],
                alignment[1] + sequenceTwo[-1]
            )
        elif index == 1:
            alignment = (
                alignment[0] + sequenceOne[-1],
                alignment[1] + "-"
            )
        elif index == 2:
            alignment = (
                alignment[0] + "-",
                alignment[1] + sequenceTwo[-1]
            )

    return score, alignment


def calculateScore(baseOne, baseTwo):
    """Returns the score for two given genetic bases"""
    if baseOne == baseTwo:
        score = SCORE_MAPPING[baseOne]
    elif baseOne == "-" or baseTwo == "-":
        score = SCORE_MAPPING["-"]
    else:
        score = SCORE_MAPPING["MISS"]

    return score


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
# Call any functions you need here, you can define them above.
# To work with the printing functions below the best alignment should be called best_alignment and its score should be called best_score. 
# The number of alignments you have checked should be stored in a variable called num_alignments.


best_score, best_alignment = alignSequences(seq1, seq2)


#-------------------------------------------------------------


# DO NOT EDIT ------------------------------------------------
# This calculates the time taken and will print out useful information 
stop = time.time()
time_taken=stop-start

# Print out the best
print('Alignments generated: '+str(num_alignments))
print('Time taken: '+str(time_taken))
print('Best (score '+str(best_score)+'):')
displayAlignment(best_alignment)

#-------------------------------------------------------------
