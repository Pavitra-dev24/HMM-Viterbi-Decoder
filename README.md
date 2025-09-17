# HMM-Viterbi-Decoder
This project parses a small discrete HMM from a text file, reads an observation sequence, runs the Viterbi dynamic programming algorithm (using log-probabilities), and prints the most likely hidden-state path and its probability.

# Files

viterbi.cpp — main source file

model.txt — example HMM model (you can edit or replace)

observations.txt — example observation sequence

# Implementation details

# Data representation

State and symbol names are stored in dynamically allocated C-style arrays of string.

Probabilities are stored in 2D arrays allocated as arrays of pointers.

Log-space arrays (log_start, log_trans, log_emit) are precomputed using log and NEG_INF for zeros.

Viterbi DP uses two 1D arrays for previous and current time-step DP scores and a 2D backpointer array of ints with dimensions T x S.

# Numerical stability

Probabilities are converted to log-space. Addition of probabilities becomes addition of log-probabilities only where appropriate (here we do max over predecessors and add emission log-probabilities).

Zero probabilities are represented by NEG_INF = -1e300. Check logic to avoid adding NEG_INF to other values directly unless guarded.
