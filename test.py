import math
import csv

# List of amino acids for matrix indexing
AMINO_ACID_SEQ = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# Factorial and combination functions
def factorial(x):
    if x == 0 or x == 1:
        return 1
    else:
        result = 1
        for i in range(2, x + 1):
            result *= i
        return result

def combination(x, y):
    if x == 0 or y == 0 or x == y:
        return 1
    dividend = factorial(x)
    divisor = factorial(x - y) * factorial(y)
    return dividend / divisor

# Print matrix in upper triangle format
def output_matrix(matrix):
    print("\t", end="")
    for amino_acid in AMINO_ACID_SEQ:
        print(amino_acid, end="\t")
    print()
    
    for i in range(len(AMINO_ACID_SEQ)):
        print(AMINO_ACID_SEQ[i], end="\t")
        for j in range(len(AMINO_ACID_SEQ)):
            if j >= i:
                print(round(matrix[i][j], 3), end="\t")
            else:
                print("-", end="\t")
        print()
    return 0

# Count occurrences of amino acid pairs in aligned sequences
def count_occurrences(aligned_sequences):
    occurrences_matrix = [[0 for _ in range(len(AMINO_ACID_SEQ))] for _ in range(len(AMINO_ACID_SEQ))]

    for i in range(len(aligned_sequences[0])):
        unique_acids_dict = {}

        # Count occurrences of each amino acid (ignoring gaps)
        for seq in aligned_sequences:
            if seq[i] == "-":
                continue
            if seq[i] not in unique_acids_dict:
                unique_acids_dict[seq[i]] = 1
            else:
                unique_acids_dict[seq[i]] += 1

        dict_keys_list = list(unique_acids_dict.keys())

        if len(unique_acids_dict) == 1:
            acid_pos = AMINO_ACID_SEQ.index(dict_keys_list[0])
            occurrences_matrix[acid_pos][acid_pos] += combination(unique_acids_dict[dict_keys_list[0]], 2)
        else:
            for i in range(len(dict_keys_list) - 1):
                for j in range(i, len(dict_keys_list)):
                    acid_i_pos = AMINO_ACID_SEQ.index(dict_keys_list[i])
                    acid_j_pos = AMINO_ACID_SEQ.index(dict_keys_list[j])

                    if i == j:
                        occurrences_matrix[acid_i_pos][acid_j_pos] += combination(unique_acids_dict[dict_keys_list[i]], 2)
                    else:
                        occurrences_matrix[acid_i_pos][acid_j_pos] += unique_acids_dict[dict_keys_list[i]] * unique_acids_dict[dict_keys_list[j]]
        
    # Apply +0.5 for all values in the upper triangle
    for i in range(len(occurrences_matrix)):
        for j in range(i + 1, len(occurrences_matrix)):
            occurrences_matrix[i][j] += 0.5
            occurrences_matrix[j][i] += 0.5

    return occurrences_matrix

# Calculate observed probability q(ij) of each amino acid pair
def calc_observed_probability(aligned_sequences, occurrences_matrix):
    row_count = len(aligned_sequences)
    col_count = len(aligned_sequences[0])
    total = 0.5 * (col_count * row_count * (row_count - 1))

    observed_probability_matrix = [[0 for _ in range(len(AMINO_ACID_SEQ))] for _ in range(len(AMINO_ACID_SEQ))]

    for i in range(len(occurrences_matrix)):
        for j in range(len(occurrences_matrix[0])):
            observed_probability = occurrences_matrix[i][j] / total
            observed_probability_matrix[i][j] = round(observed_probability, 3)

    return observed_probability_matrix

# Calculate the expected probability e(ij) of each amino acid pair
def calc_expected_probability(observed_probability_matrix):
    expected_probability_matrix = [[0 for _ in range(len(observed_probability_matrix))] for _ in range(len(observed_probability_matrix))]

    for i in range(len(observed_probability_matrix)):
        for j in range(len(observed_probability_matrix[0])):
            p_i = observed_probability_matrix[i][i]
            p_j = observed_probability_matrix[j][j]
            
            if i != j:
                expected_probability_matrix[i][j] = round(2 * p_i * p_j, 4)
            else:
                expected_probability_matrix[i][j] = round(p_i * p_i, 4)

    return expected_probability_matrix

# Calculate log-odds ratios
def calc_log_odds(observed_probability_matrix, expected_probability_matrix):
    log_odds_matrix = []
    for i in range(len(AMINO_ACID_SEQ)):
        row = []
        for j in range(len(AMINO_ACID_SEQ)):
            q_ij = observed_probability_matrix[i][j]
            e_ij = expected_probability_matrix[i][j]

            if q_ij == 0 and e_ij == 0:
                s_ij = 0
            elif q_ij == 0:
                s_ij = -10
            elif e_ij == 0:
                s_ij = 10
            else:
                s_ij = math.log2(q_ij / e_ij)

            s_ij *= 2
            s_ij = round(s_ij)

            row.append(s_ij)
        log_odds_matrix.append(row)
    return log_odds_matrix

# Save matrix to CSV
def save_matrix_to_csv(matrix, filename):
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([''] + AMINO_ACID_SEQ)

        for i in range(len(matrix)):
            row = [AMINO_ACID_SEQ[i]]
            for j in range(len(matrix)):
                if j >= i:
                    row.append(matrix[i][j])
                else:
                    row.append("")
            writer.writerow(row)
    print(f"Matrix saved to {filename}")

# Normalize sequences by padding with '-'
def normalize_sequences(sequences):
    if not sequences:
        print("Warning: No sequences provided for normalization.")
        return []
    max_length = max(len(seq) for seq in sequences)
    return [seq.ljust(max_length, '-') for seq in sequences]

# Read and parse grouped sequences from file
def read_grouped_sequences(filename):
    grouped_sequences = {}
    current_group = None
    sequences = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.endswith(':'):
                if current_group:
                    grouped_sequences[current_group] = sequences
                current_group = line
                sequences = []
            else:
                sequences.append(line)
        
        if current_group:
            grouped_sequences[current_group] = sequences

    return grouped_sequences

# Save amino acid frequencies to CSV (now calculates frequencies across all groups)
def save_amino_acid_frequencies_to_csv(all_sequences, filename):
    amino_acid_counts = {acid: 0 for acid in AMINO_ACID_SEQ}
    total_amino_acids = 0

    # Accumulate frequencies across all sequences
    for sequences in all_sequences:
        for seq in sequences:
            for acid in seq:
                if acid in AMINO_ACID_SEQ:
                    amino_acid_counts[acid] += 1
                    total_amino_acids += 1

    amino_acid_percentages = {acid: (count / total_amino_acids) * 100 for acid, count in amino_acid_counts.items()}

    # Write to CSV
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Amino Acid", "Frequency", "Percentage"])

        for acid in AMINO_ACID_SEQ:
            writer.writerow([acid, amino_acid_counts[acid], round(amino_acid_percentages[acid], 2)])

    print(f"Amino acid frequencies and percentages saved to {filename}")

# Main function
def main():
    input_file = "grouped_alignments.txt"
    sequences_by_group = read_grouped_sequences(input_file)

    # Prepare to collect all sequences for frequency calculation
    all_sequences = []

    combined_occurrences_matrix = [[0 for _ in range(len(AMINO_ACID_SEQ))] for _ in range(len(AMINO_ACID_SEQ))]

    for heading, sequences in sequences_by_group.items():
        print(f"Processing group: {heading}")

        sequences = normalize_sequences(sequences)
        occurrences_matrix = count_occurrences(sequences)

        # Collect all sequences for total frequency calculation
        all_sequences.append(sequences)

        # Accumulate occurrences matrix
        for i in range(len(combined_occurrences_matrix)):
            for j in range(len(combined_occurrences_matrix[i])):
                combined_occurrences_matrix[i][j] += occurrences_matrix[i][j]

        group_safe_heading = heading.replace(" ", "").replace("/", "").replace("(", "").replace(")", "")
        frequencies_filename = f"{group_safe_heading}_amino_acid_frequencies.csv"
        save_amino_acid_frequencies_to_csv(sequences, frequencies_filename)

    # Save frequencies for all groups combined
    save_amino_acid_frequencies_to_csv(all_sequences, 'total_amino_acid_frequencies.csv')

    save_matrix_to_csv(combined_occurrences_matrix, 'combined_occurrences_matrix.csv')

    observed_probability_matrix = calc_observed_probability(all_sequences, combined_occurrences_matrix)
    save_matrix_to_csv(observed_probability_matrix, 'combined_observed_probability_matrix.csv')

    expected_probability_matrix = calc_expected_probability(observed_probability_matrix)
    save_matrix_to_csv(expected_probability_matrix, 'combined_expected_probability_matrix.csv')

    log_odds_matrix = calc_log_odds(observed_probability_matrix, expected_probability_matrix)
    save_matrix_to_csv(log_odds_matrix, 'combined_log_odds_matrix.csv')

main()