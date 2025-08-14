import sys
import argparse
from base64 import encode
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ctypes
from fasta_reader import read_fasta

# Similar helper function as in local alignment
def create_score_matrix(alphabet: str, match_score: int = 1, mismatch_score: int = -1):
    size = len(alphabet)
    matrix = np.full((size, size), mismatch_score, dtype=np.int32)
    for i in range(size):
        matrix[i, i] = match_score
    return matrix

# Grid initialization for global alignment
def initialize_grid_global(sequence_1: str, sequence_2: str, gap_penalty: int = -1) -> np.ndarray:
    rows = len(sequence_1) + 1
    cols = len(sequence_2) + 1
    grid = np.zeros((rows, cols), dtype=np.int32)
    # Fill first column with gap penalties
    for i in range(rows):
        grid[i, 0] = i * gap_penalty
    # Fill first row with gap penalties
    for j in range(cols):
        grid[0, j] = j * gap_penalty
    return grid
    
def initialize_traceback(rows: int, cols: int) -> np.ndarray:
    return np.full((rows, cols), '', dtype='<U1')

# Python wrapper to call the C implementation of alignment
def fill_score_grid_global_c_wrapper(grid: np.ndarray, traceback: np.ndarray, sequence_1: str, sequence_2: str, score_matrix: np.ndarray, alphabet: str, gap_penalty: int):
    lib = ctypes.CDLL('./alignment_engine_global.dll') # .so on Linux/MacOS
    lib.fill_score_grid_global_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.uint8, flags='C_CONTIGUOUS'),
        ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p,
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
        ctypes.c_char_p, ctypes.c_int
    ]

    traceback_bytes = np.char.encode(traceback, 'ascii')
    seq1_bytes = sequence_1.encode('ascii')
    seq2_bytes = sequence_2.encode('ascii')
    alphabet_bytes = alphabet.encode('ascii')
    rows, cols = grid.shape

    lib.fill_score_grid_global_c(grid, traceback_bytes.view(np.uint8), rows, cols, seq1_bytes, seq2_bytes, score_matrix, alphabet_bytes, gap_penalty)

    traceback[:] = np.char.decode(traceback_bytes)
    return grid[rows-1, cols-1] # Return final alignment score

# Traceback for global alignment
def traceback_global(traceback: np.ndarray, sequence_1: str, sequence_2: str) -> tuple:
    aligned_seq1, aligned_seq2 = [], []
    i, j = traceback.shape[0] - 1, traceback.shape[1] - 1 # Start at bottom-right

    while i > 0 or j > 0:
        direction = traceback[i , j]
        if direction == 'D':
            aligned_seq1.append(sequence_1[i - 1])
            aligned_seq2.append(sequence_2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_seq1.append(sequence_1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif direction == 'L':
            aligned_seq1.append('-')
            aligned_seq2.append(sequence_2[j - 1])
            j -= 1
        else: # Reached a point with no direction (should be near top left)
            break
    return "".join(reversed(aligned_seq1)), "".join(reversed(aligned_seq2))

# Main orchestrator function
def global_alignment(sequence_1: str, sequence_2: str, alphabet: str = 'ACGT', match_score: int = 1, mismatch_score: int = -1, gap_penalty: int = -1):
    score_matrix = create_score_matrix(alphabet, match_score, mismatch_score)
    grid = initialize_grid_global(sequence_1, sequence_2, gap_penalty)
    traceback = initialize_traceback(*grid.shape)

    # Call to the C-accelerated grid filling function
    alignment_score = fill_score_grid_global_c_wrapper(grid, traceback, sequence_1, sequence_2, score_matrix, alphabet, gap_penalty)

    aligned_seq1, aligned_seq2 = traceback_global(traceback, sequence_1, sequence_2)
    
    print(f"Global Alignment Score: {alignment_score}")
    print(aligned_seq1)
    print(aligned_seq2)

    return grid
