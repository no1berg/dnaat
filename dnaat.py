import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fasta_reader import read_fasta
from global_alignment import global_alignment
from aligner import local_alignment

def main():
    parser = argparse.ArgumentParser(
        description="DNAAT: A Lightweight DNA Alignment Tool",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # --- Shared arguments for both subparsers ---
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("file", nargs='?', default=None, help="Path to the FASTA file containing two sequences to align.")
    parent_parser.add_argument("--match", type=int, default=1, help="Score for a sequence match.")
    parent_parser.add_argument("--mismatch", type=int, default=-1, help="Score for a sequence mismatch.")
    parent_parser.add_argument("--gap", type=int, default=-2, help="Gap penalty score.")
    parent_parser.add_argument("--alphabet", type=str, default="ACGTU", help="Alphabet to use for alignment.")

    # --- Global Alignment Command ---
    parser_global = subparsers.add_parser("global", parents=[parent_parser], help="Perform global alignment (Needleman-Wunsch).")
    parser_global.set_defaults(func=run_global)

    # --- Local Alignment Command ---
    parser_local = subparsers.add_parser("local", parents=[parent_parser], help="Perform local alignment (Smith-Waterman).")
    parser_local.set_defaults(func=run_local)

    args = parser.parse_args()
    args.func(args)

def get_sequences(args):
    """Reads sequences from a FASTA file or returns default examples."""
    if args.file:
        try:
            sequences = read_fasta(args.file)
            if len(sequences) < 2:
                raise ValueError("FASTA file must contain at least two sequences.")
            seq_ids = list(sequences.keys())
            seq1 = sequences[seq_ids[0]]
            seq2 = sequences[seq_ids[1]]
            print(f"Aligning '{seq_ids[0]}' and '{seq_ids[1]}' from {args.file}...")
            return seq1, seq2
        except (FileNotFoundError, ValueError) as e:
            print(f"Error: {e}")
            sys.exit(1)
    else:
        print("No FASTA file provided. Running with default example sequences.")
        if args.command == "global":
            return "GATTACA", "GCATGCU"
        else:
            return 'ACGGGGCCATACTATTATATATAATACTACGACAGTCGACTACTATATCAAAA', 'AGTTACGGGTGCCCATTTGCGCGAGCACTACGACAGTCGACTACTAGCTCAGGAAAAGG'

def run_global(args):
    """Runs the global alignment and displays the results."""
    seq1, seq2 = get_sequences(args)
    score_grid, aligned_s1, aligned_s2, score = global_alignment(seq1, seq2, alphabet=args.alphabet, match_score=args.match, mismatch_score=args.mismatch, gap_penalty=args.gap)
    
    print(f"\nGlobal Alignment Score: {score}")

    # Visualization
    score_grid_df = pd.DataFrame(score_grid, index=['-'] + list(seq1), columns=['-'] + list(seq2))
    plt.figure(figsize=(10, 8))
    plt.imshow(score_grid.astype(float), cmap='viridis', interpolation='nearest')
    plt.title("Global Alignment (Needleman-Wunsch) Score Grid")
    plt.ylabel("Sequence 1")
    plt.xlabel("Sequence 2")
    plt.colorbar(label="Alignment Score")
    plt.xticks(np.arange(len(score_grid_df.columns)), score_grid_df.columns)
    plt.yticks(np.arange(len(score_grid_df.index)), score_grid_df.index)
    
    for i in range(score_grid.shape[0]):
        for j in range(score_grid.shape[1]):
            plt.text(j, i, score_grid[i, j], ha="center", va="center", color="w", fontsize=8)
    plt.tight_layout()
    plt.show()

def run_local(args):
    """Runs the local alignment and displays the results."""
    seq1, seq2 = get_sequences(args)
    aligned_1, aligned_2, score_grid, alignment_path = local_alignment(seq1, seq2, alphabet=args.alphabet, match_score=args.match, mismatch_score=args.mismatch, gap_penalty=args.gap)

    print("\nLocal Alignment Result:")
    print(aligned_1)
    print(aligned_2)
    
    # Visualization
    score_grid_df = pd.DataFrame(score_grid, index=['-'] + list(seq1), columns=['-'] + list(seq2))
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(score_grid_df, cmap='Greys', interpolation='nearest')
    if alignment_path:
        path_rows, path_cols = zip(*alignment_path)
        ax.plot(path_cols, path_rows, color='red', linewidth=2, marker='o', markersize=3)

    ax.set_xticks(np.arange(len(score_grid_df.columns)))
    ax.set_xticklabels(score_grid_df.columns)
    ax.set_yticks(np.arange(len(score_grid_df.index)))
    ax.set_yticklabels(score_grid_df.index)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")
    fig.colorbar(im, ax=ax, label="Alignment Score")
    ax.set_title("Local Alignment Score Grid With Optimal Path")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
