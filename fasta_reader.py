# Handle the import of data from several file formats

# 1. FASTA


def read_fasta(file_path: str, validate_chars=False, validate_unique_ids=False):
    sequences = {}
    valid_chars = set('ATCGN')
    try:
        with open(file_path, 'r') as file:
            seq_id = None
            seq_lines = []
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id:
                        sequence = ''.join(seq_lines)
                        if validate_chars and not set(sequence.upper()).issubset(valid_chars):
                            raise ValueError(f"Invalid characters in sequence for ID '{seq_id}'")
                        sequences[seq_id] = sequence
                    seq_id = line[1:].strip() # Remove leading '>' and whitespace
                    if validate_unique_ids and seq_id in sequences:
                        raise ValueError(f"Duplicate sequence ID found: '{seq_id}'")
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if seq_id:
                sequence = ''.join(seq_lines)
                if valid_chars and not set(sequence.upper()).issubset(valid_chars):
                    raise ValueError(f"Invalid characters in sequence for ID '{seq_id}'")
                sequences[seq_id] = sequence

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except PermissionError:
        print(f"Error: Permission denied when accessing '{file_path}'.")
    except Exception as e:
        print(f"An error occured: {e}")
       
    return sequences