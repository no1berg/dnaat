#include <string.h> // For strlen

// A simple helper function to get the maximum of four integers
int max4(int a, int b, int c, int d) {
    int max_val = a;
    if (b > max_val) max_val = b;
    if (c > max_val) max_val = c;
    if (d > max_val) max_val = d;
    return max_val;
}

// Create a fast lookup table to get the index of a character in the alphabet
// Essentially a C-style replacement for Python's dict/hashmap
void create_alphabet_map(const char* alphabet, int* map) {
    int alphabet_size = strlen(alphabet);
    // Initialize all possible character values to -1 (not in alphabet)
    for (int i = 0; i < 256; ++i) {
        map[i] = -1;
    }
    // Map each character in the alphabet to its index
    for (int i = 0; i < alphabet_size; ++i){
        map[(unsigned char)alphabet[i]] = i;
    }
}

// Main function to perform the dynamic programing step
// All arrays (grid, traceback, score_matrix) are passed as 1D pointers
void fill_score_grid_c(
    int* grid,          // Pointer to the score grid data
    char* traceback,    // Pointer to the traceback grid data
    int rows,           // Number of rows in the grids
    int cols,           // Number of columns in the grids
    const char* seq1,   // The first sequence
    const char* seq2,   // The second sequence
    const int* score_matrix,  // Pointer to the substitution matrix data
    const char* alphabet,     // The alphabet string (e.g., "ACGT")
    int gap_penalty,          // The gap penalty
    int* max_pos              // Output array to store the [row, col] of max score
) {
    int max_score = 0;
    max_pos[0] = 0; // max_i
    max_pos[1] = 0; // max_j

    int alphabet_map[256];
    int alphabet_size = strlen(alphabet);
    create_alphabet_map(alphabet, alphabet_map);

    for (int i = 1; i < rows; ++i) {
        for (int j = 1; j < cols; ++j) {
            // Get the current characters from the sequences using the map
            int char1_idx = alphabet_map[(unsigned char)seq1[i-1]];
            int char2_idx = alphabet_map[(unsigned char)seq2[j-1]];

            // If a character is not valid, handle it here
            // Will implement at a later point

            // Calculate potential scores
            // grid[i-1, j-1] is accessed as grid[(i-1)*cols + (j-1)] in a 1D array
            int match_val = grid[(i - 1) * cols + (j - 1)] + score_matrix[char1_idx * alphabet_size + char2_idx];
            int delete_val = grid[(i - 1) * cols + j] + gap_penalty;
            int insert_val = grid[i * cols + (j - 1)] + gap_penalty;

            int score = max4(0, match_val, delete_val, insert_val);
            grid[i * cols + j] = score;

            // Update the traceback matrix based on the chosen score
            if (score == 0) {
                traceback[i * cols + j] = 'S';  // Stop
            } else if (score == match_val) {
                traceback[i * cols + j] = 'D';  // Diagonal
            } else if (score == delete_val) {
                traceback[i * cols + j] = 'U';  // Up
            } else { // score == insert_val
                traceback[i * cols + j] = 'L';  // Left
            }

            // Keep track of the highest score found so far
            if (score > max_score) {
                max_score = score;
                max_pos[0] = i;
                max_pos[1] = j;
            }
        }
    }
}