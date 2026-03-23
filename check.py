"""
This program performs pairwise sequence alignment using a substitution matrix 
and a dynamic programming approach. It takes multiple amino acid sequences 
from a FASTA file, aligns them in pairs, and combines the best-aligned sequences 
into a new one, iterating until a single sequence remains.

Input for the Program:
The input consists of two files:
1. A CSV file containing the substitution matrix, where each pair of amino acids 
   has a corresponding substitution score. 

2. A FASTA file containing multiple amino acid sequences (one per line, starting 
   with '>' followed by the sequence ID). 
   the input file must be in the following format: 
   >Alpha	
   GHAV 
   
   >Bravo	
   GHNV

Assumptions and Alerts:
1. The substitution matrix file is assumed to be named "substitution_matrix.csv".
2. The FASTA file containing the sequences is assumed to be named "files.txt".
3. The gap penalty is set to -100 for the sequence alignment.

Output format:
A final sequence that is the combination of the best-aligned sequences.

Alerts:
1. Ensure that the CSV substitution matrix file and the FASTA file are correctly formatted.


Author: TEAM Frodo 2
"""

import csv
import sys

# function to read the csv file and populate the data for subst_matrix hashmap
# example data AA: 13, AC: 8
def read_subst_matrix(filename):
    subst_matrix = {}

    # Open and read the CSV file
    with open(filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)

        # Get the list of fieldnames (amino acids)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            print("Error: CSV file does not have a header row.")
            return {}  
        fieldnames = fieldnames[1:]

        # Initialize the row index
        row_index = 0

        # Loop through each row in the CSV file
        for row in reader:
            amino1 = fieldnames[row_index]  # Get the amino acid for this row
            col_index = 0  # Initialize the column index

            # Iterate over the remaining columns (amino2 and score)
            for amino2, score in row.items():
                if col_index >= row_index + 1:  # Only read on or above the diagonal
                    if col_index == 0:  # Skip the first column as it is the row header
                        col_index += 1
                        continue

                    # Store the substitution score for the (amino1, amino2) pair
                    subst_matrix[(amino1 + amino2)] = int(score)

                col_index += 1

            row_index += 1

    # # Print out the contents of the substitution matrix
    # for key, value in subst_matrix.items():
    #     print(f"{key}: {value}")

    return subst_matrix

# function reused from previous assignment
# this function reads the file in FASTA format and store it into a hashmap with name and sequence
def read_fasta_dna_file(filename):
    amino_acids_dict = {}
    try:
        with open(filename, 'r') as file:
            name = None #declaring as a null value to check for if conditions
            actual_sequence = ""
            for line in file:
                line = line.strip()    # strips away the empty spaces in line
                if line.startswith(">"):  # This means its an header line with anme
                    if name and actual_sequence != "":  
                        amino_acids_dict[name] = actual_sequence
                        actual_sequence = ""
                    name = line[1:].strip()  # store the name after the > character
                elif line:  # if the line was not header
                    actual_sequence += line.upper()  # Accumulate the sequence
            if name and actual_sequence != "":  # Save the last sequence
                amino_acids_dict[name] = actual_sequence
    except FileNotFoundError:
        print("There is some error in the files or the file was not found")
        sys.exit()
    return amino_acids_dict

def validate_the_sequence(amino_acid_dict):
    valid_amino_acid = ['A', 'C', 'D','E','F','G','H','I','K','L','M',
                        'N','P','Q','R', 'S','T','V','W','Y', '-']  # List of valid ACID
    for sequence in amino_acid_dict.values():
        for acids in sequence:
              if acids not in valid_amino_acid:
                  return False  
    return True
    
# the dynamic programming function that will create the dp matrix and give resulant output
def create_dp_matrix(seq1, seq2, subst_matrix):
    # Initializing the DP matrix and the direction matrix
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    direction = [[''] * (n + 1) for _ in range(m + 1)]  # To store the direction (up, left, diagonal)

    # Initializing the first row and column with gap penalties (initial gap penalty is 0)
    gap_penalty = 0
    for i in range(m + 1):
        dp[i][0] = i * gap_penalty
        direction[i][0] = 'up'  # Up direction for gap in seq2
    for j in range(n + 1):
        dp[0][j] = j * gap_penalty
        direction[0][j] = 'left'  # Left direction for gap in seq1

    # Set gap penalty to -100
    gap_penalty = -100

    # Filling the rest of the DP matrix and direction matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Calculating the substitution score
            # also check if the chracter is - as found in TEST  CASE C. This was added later
            # it check if the character is - and if it is then it will give a gap penalty accordingly
            score_diagonal = dp[i - 1][j - 1] + \
            subst_matrix.get((seq1[i - 1] + seq2[j - 1]), 
                             subst_matrix.get((seq2[j - 1] + seq1[i - 1]))) if seq1[i - 1] != '-' and seq2[j - 1] != '-' else 0
            
            score_up = dp[i - 1][j] + gap_penalty
            score_left = dp[i][j - 1] + gap_penalty

            # Find the maximum score and store the corresponding direction
            dp[i][j] = max(score_diagonal, score_up, score_left)

            if dp[i][j] == score_diagonal:
                direction[i][j] = 'diagonal'  # Diagonal move means alignment
            elif dp[i][j] == score_up:
                direction[i][j] = 'up'  # Upward move means gap in seq2
            else:
                direction[i][j] = 'left'  # Leftward move means gap in seq1

    # backtracking to get the aligned sequences 
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n
    if(i == j):
        score = dp[i][j]
        row = i
        column = j
    else:
        if j > i:
            last_row_max = max(dp[i])
            row = i
            column = dp[row].index(last_row_max)
            score = dp[row][column]
        else:
            last_col_max = max([dp[i][j] for i in range(m + 1)])
            column = j
            row = 0
            for r in range(m + 1):
                if dp[r][column] == last_col_max:
                    row = r
            score = dp[row][column]
            
    while row > 0 or column > 0:
        if direction[row][column] == 'diagonal':
            # Diagonal move: both characters are aligned
            aligned_seq1.append(seq1[row - 1])
            aligned_seq2.append(seq2[column - 1])
            row -= 1
            i -= 1
            column -= 1
            j -= 1
        elif direction[row][column] == 'up':
            # Up move: seq2 has a gap
            aligned_seq1.append(seq1[row - 1])
            aligned_seq2.append('-')
            row -= 1
            i -= 1

        elif direction[row][column] == 'left':
            # Left move: seq1 has a gap
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[column - 1])
            column -= 1
            j -= 1

    # Reverse the aligned sequences since we traced them backwards 
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    # If there are any remaining characters in seq1 or seq2
    while abs(i - row) > 0:
        
        aligned_seq1.append(seq1[-i])
        aligned_seq2.append('-')
        i -= 1

    while abs(j - column) > 0:
        aligned_seq1.append('-')
        aligned_seq2.append(seq2[-j])
        j -= 1

    
    return dp, score, ''.join(aligned_seq1), ''.join(aligned_seq2), direction

# Function to combine two aligned sequences into one
# if there is a mistach like GHAV and GHNV 
# this function will consult the subst_matrix and choose the candiate who have 
# higher value for self match, as it is more likely that other postison will also be that
def mix_sequences(seq1, seq2, subst_matrix):
    combined = []
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            combined.append(s1)  # Match
        elif s1 == '-':
            combined.append(s2)  # Insert from seq2
        elif s2 == '-':
            combined.append(s1)  # Insert from seq1
        else:
            # If characters mismatch, consult the substitution matrix
            score1 = subst_matrix.get((s1 + s1), 0)  
            score2 = subst_matrix.get((s2 + s2), 0)  

            if score1 >= score2:
                combined.append(s1)
            else:
                combined.append(s2)
    return ''.join(combined)

def main():
    # Reading the data from csv
    subst_matrix = read_subst_matrix("combined_log_odds_matrix.csv")

    # Reading sequences from the FASTA file
    dict_sequences = read_fasta_dna_file("sequence.txt")
    bool_decision = validate_the_sequence(dict_sequences)
    if not bool_decision:
        print("Your Sequences are invalid check your file again")
        sys.exit()

    # Convert dict_sequences into a list of tuples (id, sequence) for easy iteration
    sequences = list(dict_sequences.items())

    # List to store results of alignments
    alignments = []
    all_aligned_values = {}

    while len(sequences) > 1:
        highest_score = float('-inf')
        best_alignment = None
        best_indices = None
        
        # Iterate through all pairs of sequences
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                id1, seq1 = sequences[i]
                id2, seq2 = sequences[j]

                # Perform sequence alignment
                matrix, score, seq1a, seq2a, direction = create_dp_matrix(seq1, seq2, subst_matrix)

                # Check if this is the best score
                if score > highest_score:
                    highest_score = score
                    best_alignment = (id1, seq1a, id2, seq2a, score)
                    best_indices = (i, j)

        # If no alignment was found, break the loop
        if best_alignment is None:
            break
        
        # print("Round Best Score:", highest_score)
        # print("Sequence 1 ID:", best_alignment[0])
        # print("Aligned Sequence 1:", best_alignment[1])
        # print("Sequence 2 ID:", best_alignment[2])
        # print("Aligned Sequence 2:", best_alignment[3])

        # Combining the two aligned sequences
        combined_seq = mix_sequences(best_alignment[1], best_alignment[3], subst_matrix)
        combined_id = f"{best_alignment[0]}_{best_alignment[2]}"

        # Remove the aligned sequences from the hashmap and adding the new one
        if best_indices is not None:
            i, j = best_indices
            del sequences[j]  
            del sequences[i]

        # Adding the combined sequence back to the list
        sequences.append((combined_id, combined_seq))

        # Storing the alignment for reporting purpose
        if best_alignment[1] != combined_seq:
            alignments.append((best_alignment[1]))
        else:
            alignments.append(best_alignment[3])
        
        


    for key, values in sequences:
        alignments.append(values) 
        # adding the last combined sequence to the alignments
        print("\nFinal Signal Alignment Found:")
        print("--------------\n")
        print(f"Combined Fasta Label: {key}")
        print("--------------\n")
        print(f"Combined  Sequence: {values}")
        print("--------------\n")
    
    print("Alignment of the Sequences")
    print("--------------\n")
    for values in alignments:
        print(values)
        

if __name__ == "__main__":
    main()