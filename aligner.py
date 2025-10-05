import argparse
import sys
import ctypes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fasta_reader import read_fasta
from typing import Dict, Any, Tuple

def create_score_matrix(alphabet: str, match_score: int = 1, mismatch_score: int = -1) -> Tuple[np.ndarray, Dict[str, int]]:
    """Creates a score matrix for the given alphabet."""
    size = len(alphabet)
    matrix = np.full((size, size), mismatch_score, dtype=int)
    for i in range(size):
        matrix[i, i] = match_score
    return matrix, {char: idx for idx, char in enumerate(alphabet)}

def initialize_grid(sequence_1: str, sequence_2: str) -> np.ndarray:
    """Initializes the grid for local alignment."""
    rows = len(sequence_1) + 1
    cols = len(sequence_2) + 1
    grid = np.zeros((rows, cols), dtype=int) # Initialize with all zeroes, allowing start anywhere
    return grid

def initialize_traceback(rows: int, cols: int) -> np.ndarray:
    """Initializes the traceback matrix."""
    return np.full((rows, cols), '', dtype='<U1') # U1 = Unicode string of length 1

def fill_score_grid_c_wrapper(grid: np.ndarray, traceback: np.ndarray, sequence_1: str, sequence_2: str, score_matrix: np.ndarray, alphabet: str, gap_penalty: int) -> tuple:
    """Python wrapper to call the C implementation of fill_score_grid."""
    # Load the shared library, .so (Linux/MacOS) and .dll (Windows)
    lib = ctypes.CDLL('./alignment_engine.dll')

    # Define the function signatures (argument types)
    lib.fill_score_grid_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.uint8, flags='C_CONTIGUOUS'), # For char array
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_char_p,
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS'),
        ctypes.c_char_p,
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS')
    ]

    # The traceback matrix in Python was '<U1', but here it must be bytes for C
    traceback_bytes = np.char.encode(traceback, 'ascii')

    # Prepare the sequences for C (encode to bytes)
    seq1_bytes = sequence_1.encode('ascii')
    seq2_bytes = sequence_2.encode('ascii')
    alphabet_bytes = alphabet.encode('ascii')

    # Prepare output array for max_pos
    max_pos_out = np.zeros(2, dtype=np.int32)
    rows, cols = grid.shape

    # Ensure arrays are C-contiguous and have the correct d-type for the C function
    grid32 = np.ascontiguousarray(grid, dtype=np.int32)
    score_matrix32 = np.ascontiguousarray(score_matrix, dtype=np.int32)

    # Call the C function
    lib.fill_score_grid_c(
        grid32,
        traceback_bytes.view(np.uint8),
        rows,
        cols,
        seq1_bytes,
        seq2_bytes,
        score_matrix32,
        alphabet_bytes,
        gap_penalty,
        max_pos_out
    )

    # Update the original arrays with results from C
    grid[:] = grid32[:]
    traceback[:] = np.char.decode(traceback_bytes)

    # Return the max position as a tuple
    return tuple(max_pos_out)

def traceback_local(traceback: np.ndarray, sequence_1: str, sequence_2: str, start_pos: tuple, grid: np.ndarray) -> Tuple[str, str, int, list]:
    """Performs traceback for local alignment."""
    aligned_sequence1 = []
    aligned_sequence2 = []
    i, j = start_pos
    path = [(i, j)]
    score = grid[i, j]

    while traceback[i, j] != '' and traceback[i, j] != 'S':
        direction = traceback[i ,j]
        if direction == 'D':
            aligned_sequence1.append(sequence_1[i-1])
            aligned_sequence2.append(sequence_2[j-1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_sequence1.append(sequence_1[i-1])
            aligned_sequence2.append('-')
            i -= 1
        elif direction == 'L':
            aligned_sequence1.append('-')
            aligned_sequence2.append(sequence_2[j-1])
            j -= 1
        else:
            break
        path.append((i, j))

    return ''.join(reversed(aligned_sequence1)), ''.join(reversed(aligned_sequence2)), score, path

def local_alignment(sequence_1: str, sequence_2: str, alphabet: str = 'ACGT', match_score: int = 1, mismatch_score: int = -1, gap_penalty: int = -1) -> Dict[str, Any]:
    """Performs local alignment on two sequences."""
    score_matrix, index_map = create_score_matrix(alphabet, match_score, mismatch_score)
    grid = initialize_grid(sequence_1, sequence_2)
    traceback = initialize_traceback(*grid.shape)
    max_pos = fill_score_grid_c_wrapper(grid, traceback, sequence_1, sequence_2, score_matrix, alphabet, gap_penalty)
    aligned_sequence1, aligned_sequence2, score, path = traceback_local(traceback, sequence_1, sequence_2, max_pos, grid)

    return {
        "score_grid": grid,
        "aligned_s1": aligned_sequence1,
        "aligned_s2": aligned_sequence2,
        "score": score,
        "alignment_path": path
    }