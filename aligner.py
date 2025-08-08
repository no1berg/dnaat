import argparse
import sys
import ctypes
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fasta_reader import read_fasta

def create_score_matrix(alphabet: str, match_score: int = 1, mismatch_score: int = -1):
    size = len(alphabet)
    matrix = np.full((size, size), mismatch_score, dtype=int)
    for i in range(size):
        matrix[i, i] = match_score
    return matrix, {char: idx for idx, char in enumerate(alphabet)}

def initialize_grid(sequence_1: str, sequence_2: str, gap_penalty: int = -1) -> np.ndarray:
    rows = len(sequence_1) + 1
    cols = len(sequence_2) + 1
    grid = np.zeros((rows, cols), dtype=int) # Initialize with all zeroes, allowing start anywhere
    return grid


def initialize_traceback(rows: int, cols: int) -> np.ndarray:
    return np.full((rows, cols), '', dtype='<U1') # U1 = Unicode string of length 1

def fill_score_grid_c_wrapper(grid: np.ndarray, traceback: np.ndarray, sequence_1: str, sequence_2: str, score_matrix: np.ndarray, alphabet: str, gap_penalty: int):
    "Python wrapper to call the C implementation of fill_score_grid."
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


def traceback_local(traceback: np.ndarray, sequence_1: str, sequence_2: str, start_pos: tuple) -> tuple:
    aligned_sequence1 = []
    aligned_sequence2 = []
    i, j = start_pos
    path = [(i, j)]

    while traceback[i, j] != '':
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
        path.append((i, j))

    return ''.join(reversed(aligned_sequence1)), ''.join(reversed(aligned_sequence2)), path
   
def local_alignment(sequence_1: str, sequence_2: str, alphabet: str = 'ACGT', match_score: int = 1, mismatch_score: int = -1, gap_penalty: int = -1):
    score_matrix, index_map = create_score_matrix(alphabet, match_score, mismatch_score)
    grid = initialize_grid(sequence_1, sequence_2, gap_penalty)
    traceback = initialize_traceback(*grid.shape)
    max_pos = fill_score_grid_c_wrapper(grid, traceback, sequence_1, sequence_2, score_matrix, alphabet, gap_penalty)
    aligned_sequence1, aligned_sequence2, path = traceback_local(traceback, sequence_1, sequence_2, max_pos)
    return aligned_sequence1, aligned_sequence2, grid, path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform local DNA sequence alignment.")
    parser.add_argument("file", nargs='?', default=None, help="Path to the FASTA file containing two sequences to align.")
    
    args = parser.parse_args()

    if args.file:
        try:
            sequences = read_fasta(args.file)
            if len(sequences) < 2:
                raise ValueError("FASTA file must contain at least two sequences.")
            
            # Extract the first two sequences from the file
            seq_ids = list(sequences.keys())
            sequence_1 = sequences[seq_ids[0]]
            sequence_2 = sequences[seq_ids[1]]
            print(f"Aligning '{seq_ids[0]}' and '{seq_ids[1]}' from {args.file}...")

        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        # Default example sequences if no file is provided
        print("No FASTA file provided. Running with default example sequences.")
        sequence_1 = 'ACGGGGCCATACTATTATATATAATACTACGACAGTCGACTACTATATCAAAA'
        sequence_2 = 'AGTTACGGGTGCCCATTTGCGCGAGCACTACGACAGTCGACTACTAGCTCAGGAAAAGG'

    # Perform local alignment
    aligned_1, aligned_2, score_grid, alignment_path, s1, s2 = local_alignment(sequence_1, sequence_2)
    
    print("\nLocal Alignment Result:")
    print(aligned_1)
    print(aligned_2)

    # Create DataFrame for visualization
    score_grid_df = pd.DataFrame(score_grid)

    # Plot 1. Full score grid heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(score_grid_df, cmap='Greys')

    # Overlay alignment path
    path_rows, path_cols = zip(*alignment_path)
    ax.plot(path_cols, path_rows, color='red', linewidth=2, marker='o', markersize=3)

    # Set axis labels
    labels_rows = [''] + list(s1)
    labels_cols = [''] + list(s2)

    ax.set_yticks(np.arange(len(labels_rows)))
    ax.set_yticklabels(labels_rows)

    ax.set_xticks(np.arange(len(labels_cols)))
    ax.set_xticklabels(labels_cols)
    
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")

    fig.colorbar(im, ax=ax)
    ax.set_title("Local Alignment Score Grid With Optimal Path")

    plt.tight_layout()
    plt.show()