#include <string.h> // For strlen

// Helper function to get the maximum of three integers
int max3(int a, int b, int c) {
    int max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    return max_val;
}

// Use the same fast lookup table as the local alignment engine
void create_alphabet_map(const char* alphabet, int* map) {
    int alphabet_size = strlen(alphabet);
    for (int i = 0; i < 256; ++i) {
        map[i] = -1;
    }
    for (int i = 0; i < alphabet_size; ++i) {
        map[(unsigned char)alphabet[i]] = i;
    }
}

// C function for the Needleman-Wunsch scoring grid
void fill_score_grid_global_c(
    int* grid,               // Pointer to the score grid data
    char* traceback,         // Pointer to the traceback grid data
    int rows,                // Number of rows in the grids
    int cols,                // Number of columns in the grids
    const char* seq1,        // The first sequence
    const char* seq2,        // The second sequence
    const int* score_matrix,    // Pointer to the substitution matrix
    const char* alphabet,       // The alphabet string
    int gap_penalty             // The gap penalty
) {
    int alphabet_map[256];
    int alphabet_size = strlen(alphabet);
    create_alphabet_map(alphabet, alphabet_map);

    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            int char1_idx = alphabet_map[(unsigned char)seq1[i - 1]];
            int char2_idx = alphabet_map[(unsigned char)seq2[j - 1]];

            int match_val = grid[(i - 1) * cols + (j - 1)] + score_matrix[char1_idx * alphabet_size + char2_idx];
            int delete_val = grid[(i - 1) * cols + j] + gap_penalty;
            int insert_val = grid[i * cols + (j - 1)] + gap_penalty;

            int score = max3(match_val, delete_val, insert_val);
            grid[i * cols + j] = score;

            // Update traceback
            // Simple tie-breaking: Diagonal > Up > Left
            if (score == match_val) {
                traceback[i * cols + j] = 'D';  // Diagonal
            } else if (score == delete_val) {
                traceback[i * cols + j] = 'U';  // Up
            } else {
                traceback[i * cols + j] = 'L';  // Left
            }
        }
    }
}