import timeit
import numpy as np
import ctypes
import random
import string

# ===================================================================
# 1. PASTE YOUR ORIGINAL PYTHON & C-WRAPPER FUNCTIONS HERE
# ===================================================================

# --- Original Python Function ---
def create_score_matrix(alphabet: str, match_score: int = 1, mismatch_score: int = -1):
    size = len(alphabet)
    matrix = np.full((size, size), mismatch_score, dtype=np.int32)
    for i in range(size):
        matrix[i, i] = match_score
    return matrix, {char: idx for idx, char in enumerate(alphabet)}

def fill_score_grid_python(grid: np.ndarray, traceback: np.ndarray, sequence_1: str, sequence_2: str, score_matrix: np.ndarray, index_map: dict, gap_penalty: int):
    max_score = 0
    max_pos = (0, 0)
    for i in range(1, len(sequence_1) + 1):
        for j in range(1, len(sequence_2) + 1):
            match = grid[i-1, j-1] + score_matrix[index_map[sequence_1[i-1]], index_map[sequence_2[j-1]]]
            delete = grid[i-1, j] + gap_penalty
            insert = grid[i, j-1] + gap_penalty
            score = max(0, match, delete, insert)
            grid[i, j] = score
            if score > max_score:
                max_score = score
                max_pos = (i, j)
    return max_pos

# --- C Wrapper Function ---
def fill_score_grid_c_wrapper(grid: np.ndarray, traceback: np.ndarray, sequence_1: str, sequence_2: str, score_matrix: np.ndarray, alphabet: str, gap_penalty: int):
    try:
        lib = ctypes.CDLL('./alignment_engine.dll')
    except OSError:
        print("Error: Could not find the compiled C library (alignment_engine.so/dll).")
        print("Please make sure you have compiled the C code first.")
        return None

    lib.fill_score_grid_c.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS"),
        ctypes.c_int, ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),
        ctypes.c_char_p, ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
    ]
    
    traceback_bytes = np.char.encode(traceback, 'ascii')
    seq1_bytes = sequence_1.encode('ascii')
    seq2_bytes = sequence_2.encode('ascii')
    alphabet_bytes = alphabet.encode('ascii')
    max_pos_out = np.zeros(2, dtype=np.int32)
    rows, cols = grid.shape
    grid32 = np.ascontiguousarray(grid, dtype=np.int32)
    score_matrix32 = np.ascontiguousarray(score_matrix, dtype=np.int32)

    lib.fill_score_grid_c(grid32, traceback_bytes.view(np.uint8), rows, cols, seq1_bytes, seq2_bytes, score_matrix32, alphabet_bytes, gap_penalty, max_pos_out)
    return tuple(max_pos_out)

# ===================================================================
# 2. TEST SETUP AND EXECUTION
# ===================================================================
if __name__ == "__main__":
    print("🔬 Setting up performance test...")

    # --- Test Parameters ---
    SEQ_LENGTH = 500 # Use a longer sequence to see a meaningful difference
    ALPHABET = 'ACGT'
    GAP_PENALTY = -1
    MATCH_SCORE = 1
    MISMATCH_SCORE = -1
    NUM_RUNS = 10 # Number of times to run each function for timing

    # --- Generate long random sequences for a realistic test ---
    sequence_1 = ''.join(random.choices(ALPHABET, k=SEQ_LENGTH))
    sequence_2 = ''.join(random.choices(ALPHABET, k=SEQ_LENGTH))
    print(f"Generated two random sequences of length {SEQ_LENGTH}.")

    # --- Prepare common inputs outside the timer ---
    score_matrix, index_map = create_score_matrix(ALPHABET, MATCH_SCORE, MISMATCH_SCORE)
    rows = len(sequence_1) + 1
    cols = len(sequence_2) + 1
    
    # We create fresh grids for each run to avoid timing the reset operation
    
    print(f"Benchmarking over {NUM_RUNS} runs...")
    print("-" * 30)

    # --- Time the Python version ---
    python_setup = f"""
from __main__ import fill_score_grid_python, sequence_1, sequence_2, score_matrix, index_map, GAP_PENALTY, rows, cols
import numpy as np
# Create fresh inputs for each timed run
grid = np.zeros((rows, cols), dtype=np.int32)
traceback = np.full((rows, cols), '', dtype='<U1')
"""
    python_stmt = "fill_score_grid_python(grid, traceback, sequence_1, sequence_2, score_matrix, index_map, GAP_PENALTY)"
    python_time = timeit.timeit(stmt=python_stmt, setup=python_setup, number=NUM_RUNS)
    print(f"🐍 Python version took: {python_time:.4f} seconds")

    # --- Time the C version ---
    c_setup = f"""
from __main__ import fill_score_grid_c_wrapper, sequence_1, sequence_2, score_matrix, ALPHABET, GAP_PENALTY, rows, cols
import numpy as np
# Create fresh inputs for each timed run
grid = np.zeros((rows, cols), dtype=np.int32)
traceback = np.full((rows, cols), '', dtype='<U1')
"""
    c_stmt = "fill_score_grid_c_wrapper(grid, traceback, sequence_1, sequence_2, score_matrix, ALPHABET, GAP_PENALTY)"
    c_time = timeit.timeit(stmt=c_stmt, setup=c_setup, number=NUM_RUNS)
    print(f"⚙️ C version took:      {c_time:.4f} seconds")
    
    print("-" * 30)

    # --- Calculate and report the results ---
    if c_time > 0:
        speedup = python_time / c_time
        print(f"✅ Performance Boost: The C version is {speedup:.2f}x faster!")
    else:
        print("C version ran too fast to measure a significant speedup.")
