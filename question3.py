#!/usr/bin/python
import numpy as np


def NJ(fileName):
    """ Merges a species Distance Matrix using the Neighbour Joining algorithm."""
    species, distanceMatrix = readData(fileName)
    speciesNum = len(species)

    # speciesNum also corresponds to the dimension of the matrix
    while speciesNum >= 2:
        rowSums = [sum(row) for row in distanceMatrix]

        print("\nDISTANCE + (row scores):")
        printMatrix(species, distanceMatrix, rowSums)

        # Here we generate the Q scores for the current distance matrix to
        # find the two species with the lowest Q score.
        qScoreMatrix = generateQScores(speciesNum, rowSums, distanceMatrix)
        lowestSpecies = findLowestScore(qScoreMatrix)

        print("\nQ SCORES:")
        printMatrix(species, qScoreMatrix)

        # We merge the two lowest species and create a new
        # corresponding distance matrix.
        species, distanceMatrix = mergeSpecies(speciesNum, species, lowestSpecies, distanceMatrix)
        speciesNum = len(species)


# region Input/Output
def readData(fileName):
    """Extracts species and their distance matrix from file given."""
    with open(fileName, "r") as matrixFile:
        lines = [line.strip("\n").split(" ") for line in matrixFile]

    species = lines.pop(0)[1:]
    distanceMatrix = np.array(
        [
            line[1:]
            for line in lines
        ],
        dtype=float
    )
    return species, distanceMatrix


def printMatrix(side, matrix, rightAppend=None):
    """Prints the given matrix, decorating with the values given."""
    if rightAppend is None:
        rightAppend = ["" for i in side]
    else:
        rightAppend = [f"({i})" for i in rightAppend]

    # This converts the default string representation of a numpy array
    # into one that is slightly cleaner and stores each line
    # separately in a list.
    matrixStrings = [row[2:].strip("]") for row in str(matrix).split("\n")]

    # Here we calculate the spacers necessary to place thing vaguely in
    # alignment with each other.
    topSpacer = int(len(matrixStrings[0])/len(side))*" "
    hSpacer = (int(max([len(sideValue) for sideValue in side])) + 2)

    # We generate the individual rows by combing the sides, spacers
    # and matrix rows into indivual row strings.
    topRow = hSpacer*" " + "| " + "".join([sideValue + topSpacer for sideValue in side])
    rows = [topRow] + [
        side[i] + (hSpacer-len(side[i]))*" " + "| " + matrixString + hSpacer*" " + str(rightAppend[i])
        for i, matrixString in enumerate(matrixStrings)
    ]
    # We then print these joined by newlines.
    print("\n".join(rows))

# endregion Input/Output


# region Q Scores
def generateQScores(speciesNum, rowSums, distanceMatrix):
    """Creates the QScore Matrix for a given species Matrix,"""
    qScoreMatrix = np.empty(distanceMatrix.shape, dtype=float)

    for i in range(speciesNum):
        for j in range(speciesNum):
            # Check if we're calculating a score for one species.
            if j == i:
                qScore = 0
            # Check if we've calculated the score before due to it
            # being a symmetrical matrix
            elif j < i:
                qScore = qScoreMatrix[j][i]
            # If not calculate the Q score.
            else:
                qScore = calculateQScore(
                    (i, j),
                    speciesNum,
                    rowSums,
                    distanceMatrix
                )

            qScoreMatrix[i][j] = qScore

    return qScoreMatrix


def calculateQScore(speciesIndexes, speciesNum, rowSums, distanceMatrix):
    """Calculate the Q score for two given species."""
    # Q(a,b) = (r-1)*d(a, b) - sum(d(a, i)) - sum(d(b, i))
    qScore = (
        (speciesNum - 1)*distanceMatrix[speciesIndexes]
        - sum([
            rowSums[speciesIndex]
            for speciesIndex in speciesIndexes
        ])
    )

    return str(qScore)


def findLowestScore(qMatrix):
    """Returns the species with the lowest score"""
    lowestSpecies = np.unravel_index(
        np.argmin(np.triu(qMatrix), axis=None),
        qMatrix.shape
    )

    # We return the sorted list to make understanding
    # referencing easier later
    return tuple(sorted(lowestSpecies))

# endregion Q Scores


# region Merge Species
def mergeSpecies(speciesNum, species, lowestSpecies, distanceMatrix):
    """Merges two given species to return a reduced distance matrix."""
    for i in range(speciesNum):
        if i not in lowestSpecies:
            newDistance = calculateNewDistance(i, lowestSpecies, distanceMatrix)

            # The lowest indexed lowest specie will become the merged specie so
            # insert the new distance for this merged specie.
            distanceMatrix[(i, lowestSpecies[0])] = newDistance
            distanceMatrix[(lowestSpecies[0], i)] = newDistance

    # Delete all entries for the highest indexed lowest specie.
    for i in range(2):
        distanceMatrix = np.delete(distanceMatrix, (lowestSpecies[1]), axis=i)

    # Generate the merged species' name and insert in the correct place
    mergedSpecies = "".join([species[lowestSpecie] for lowestSpecie in lowestSpecies])
    species[lowestSpecies[0]] = mergedSpecies
    # Remove the highest indexed lowest specie
    species.pop(lowestSpecies[1])

    return species, distanceMatrix


def calculateNewDistance(otherSpecie, lowestSpecies, distanceMatrix):
    """Calculates the distance between the merged species and the other species."""
    # d(ab, c) = (d(a,i) + d(b,j) - d(a,b))/2
    distanceScore = (
        sum([
            distanceMatrix[(otherSpecie, lowestSpecie)]
            for lowestSpecie in lowestSpecies
        ])
        - distanceMatrix[lowestSpecies]
    )/2

    return distanceScore

# endregion Merge Species

