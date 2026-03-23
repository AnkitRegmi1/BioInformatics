''' Psuedocode:
We start by importing the csv library.
Inside main() we do the following steps:
We set the file paths and gap penalty.
We then call the read_substitution() to read the substitution matrix.
We call read_fasta() to read the sequences from the FASTA file.
we call progressive_align() with the sequences, substitution matrix, and gap penalty.

Inside progressive_align() we perform the following steps:
We take the first sequence as the reference.

For each subsequent sequence we perform the following steps:
global_alignment() is called to align it with the reference.
The reference is updated with the new alignment.
All previously aligned sequences are adjusted to the new length.
Finally, adjust_to_significant_length() is called to ensure all sequences are the same length.

The code then goes back in main() and calls format_alignment() to format the aligned sequences.



'''

import csv


def read_matrix(file_path):
    
    substitution_matrix = {}
    with open(file_path, "r") as file:
        reader = csv.reader(file)
        amino_acids = next(reader)[1:]  # Readign the header row
        for row in reader:
            if not row or len(row) < len(amino_acids) + 1:
                continue
            row_header = row[0]
            scores = row[1:]
            for col_header, score in zip(amino_acids, scores):
                substitution_matrix[(row_header, col_header)] = int(score) if score.strip() else 0
    return substitution_matrix


def read_fasta(file_path):
    
    sequences = {}
    with open(file_path, "r") as file:
        seq_id = None
        seq = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id.split()[0]] = seq  
                seq_id = line[1:]
                seq = ""
            else:
                seq += line
        if seq_id:
            sequences[seq_id.split()[0]] = seq  
    return sequences


def global_alignment(seq1, seq2, matrix, gap_penalty):
    
    m, n = len(seq1), len(seq2)
    # Initializing the DP table (we could optimize memory here, but keeping it simple for now)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[''] * (n + 1) for _ in range(m + 1)]

    # Initializing the DP matrix
    for i in range(1, m + 1):
        dp[i][0] = i * gap_penalty
        traceback[i][0] = "up"
    for j in range(1, n + 1):
        dp[0][j] = j * gap_penalty
        traceback[0][j] = "left"

 
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + matrix.get(
                (seq1[i - 1], seq2[j - 1]), gap_penalty
            )
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(match, delete, insert)

            if dp[i][j] == match:
                traceback[i][j] = "diagonal"
            elif dp[i][j] == delete:
                traceback[i][j] = "up"
            else:
                traceback[i][j] = "left"

    # Traceback
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == "diagonal":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == "up":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()
    return dp[m][n], ''.join(aligned_seq1), ''.join(aligned_seq2)


def adjust_length(aligned_sequences):
    
    max_length = max(len(seq.rstrip("-")) for seq in aligned_sequences.values())

    # Adjusting all sequences to the max significant length so that all the sequences are of the same length
    for seq_id in aligned_sequences:
        aligned_sequences[seq_id] = aligned_sequences[seq_id][:max_length].ljust(max_length, "-")

    return aligned_sequences


def progressive_align(sequences, substitution_matrix, gap_penalty):
   
    aligned_sequences = {seq_id: seq for seq_id, seq in sequences.items()}

    # Starting with the first sequence as the reference
    reference_seq= list(aligned_sequences.keys())[0]
    ref_seq = aligned_sequences[reference_seq]

    for seq_id in list(aligned_sequences.keys())[1:]:
        seq = aligned_sequences[seq_id]
        _, aligned_ref, aligned_seq = global_alignment(ref_seq, seq, substitution_matrix, gap_penalty)
        ref_seq = aligned_ref

        # Ensuring that all sequences are adjusted to the new reference alignment lengt
        max_length = max(len(ref_seq), len(aligned_seq))
        for existing_id in aligned_sequences.keys():
            aligned_sequences[existing_id] = aligned_sequences[existing_id].ljust(max_length, "-")
        aligned_sequences[seq_id] = aligned_seq.ljust(max_length, "-")

    
    aligned_sequences = adjust_length(aligned_sequences)

    return aligned_sequences


def format_alignmnt(aligned_sequences, line_width=60):
    
    maximum_length = max(len(seq_id) for seq_id in aligned_sequences)
    formatted_output = ""

    # Breaking alignmnt into chunks of line_width for each sequence
    num_blocks = (len(next(iter(aligned_sequences.values()))) + line_width - 1) // line_width
    for block_idx in range(num_blocks):
        start = block_idx * line_width
        end = start + line_width

        # Adding each sequence's chunk
        for seq_id, aligned_seq in aligned_sequences.items():
            seq_chunk = aligned_seq[start:end]
            formatted_output += f"{seq_id.ljust(maximum_length)}  {seq_chunk}\n"

        formatted_output += "\n"  

    return formatted_output


def main():
    substitution_matrix_file = "combined_log_odds_matrix.csv"
    fasta_file_path = "sequence.txt"
    # Not sure if the gap penalty should be -50, but we sticked with this as it is larger than any score in our log_oods_matrix(which is our substitution matrix)
    gap_penalty = -50

    # Reading the substitution matrix and sequences
    substitution_matrix = read_matrix(substitution_matrix_file)
    sequences = read_fasta(fasta_file_path)

    # Performing progressive alignment
    aligned_sequences = progressive_align(sequences, substitution_matrix, gap_penalty)

    # Format and output the alignmnt in phylogeny.fr style so that we can compare for the similarity between the original alignment and the customised alignment
    formatted_alignment = format_alignmnt(aligned_sequences)
    print("\nMultiple Sequence Alignment:\n")
    print(formatted_alignment)


if __name__ == "__main__":
    main()
