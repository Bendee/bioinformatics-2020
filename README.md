# Bioinformatics Coursework 2020
(2020) Bioinformatics submodule coursework submission for the Computational Thinking module at Durham University

## Description
This coursework consisted of three parts.
The first two questions involved implementing an algorithm for aligning DNA sequences, and the third required the implementation of an algorithm that would allow building [phylogenetic trees](https://wikipedia.org/wiki/Phylogenetic_tree).

## Questions
1. Computes the optimal alignment of two DNA sequences using recursion.
2. Computes optimal local alignment of two sequences using a backtracking matrix.
3. Implements the [Neighbour Joining algorithm](https://wikipedia.org/wiki/Neighbor_joining) to continuously cluster species together.

### Alignment Scoring
The alignment scoring used for Q1 and Q2 is as follows:
Alignment | Score
----------|------
Matching 'A' bases | +3
Matching 'C' bases | +2
Matching 'G' bases | +1
Matching 'T' bases | +2
Mismatching bases | -3
Gap | -4
