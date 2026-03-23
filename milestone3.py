#Milestone 3 Question 3

#We used the grantham_matrix
#function to look up the distance bw two amino acids in the matrix
#analyze alignmnet, which will check for identical, conserved and variable
#alignments will be processed by analyzing data from both original and customized alignment txt files
#results will be calculated, bar chart to comapre the result

import matplotlib.pyplot as plt
#library for creating visualizations

# creates the grantham distance matrix for amino acid comparison
def grantham_matrix():
    grantham = {
        ('R', 'A'): 112, ('R', 'N'): 86, ('R', 'D'): 96, ('R', 'C'): 180, ('R', 'Q'): 43,
        ('R', 'E'): 54, ('R', 'G'): 125, ('R', 'H'): 29, ('R', 'I'): 97, ('R', 'L'): 102,
        ('R', 'K'): 26, ('R', 'M'): 91, ('R', 'F'): 97, ('R', 'P'): 103, ('R', 'S'): 110,
        ('R', 'T'): 71, ('R', 'W'): 101, ('R', 'Y'): 77, ('R', 'V'): 96, ('N', 'A'): 111,
        ('N', 'D'): 23, ('N', 'C'): 139, ('N', 'Q'): 46, ('N', 'E'): 42, ('N', 'G'): 80,
        ('N', 'H'): 68, ('N', 'I'): 149, ('N', 'L'): 153, ('N', 'K'): 94, ('N', 'M'): 142,
        ('N', 'F'): 158, ('N', 'P'): 91, ('N', 'S'): 46, ('N', 'T'): 65, ('N', 'W'): 174,
        ('N', 'Y'): 143, ('N', 'V'): 133, ('D', 'A'): 126, ('D', 'C'): 154, ('D', 'Q'): 61,
        ('D', 'E'): 45, ('D', 'G'): 94, ('D', 'H'): 81, ('D', 'I'): 168, ('D', 'L'): 172,
        ('D', 'K'): 101, ('D', 'M'): 160, ('D', 'F'): 177, ('D', 'P'): 108, ('D', 'S'): 65,
        ('D', 'T'): 85, ('D', 'W'): 181, ('D', 'Y'): 160, ('D', 'V'): 152, ('C', 'A'): 195,
        ('C', 'Q'): 154, ('C', 'E'): 170, ('C', 'G'): 159, ('C', 'H'): 174, ('C', 'I'): 198,
        ('C', 'L'): 198, ('C', 'K'): 202, ('C', 'M'): 196, ('C', 'F'): 205, ('C', 'P'): 169,
        ('C', 'S'): 112, ('C', 'T'): 149, ('C', 'W'): 215, ('C', 'Y'): 194, ('C', 'V'): 192,
        ('Q', 'A'): 91, ('Q', 'E'): 29, ('Q', 'G'): 87, ('Q', 'H'): 24, ('Q', 'I'): 109,
        ('Q', 'L'): 113, ('Q', 'K'): 53, ('Q', 'M'): 101, ('Q', 'F'): 116, ('Q', 'P'): 76,
        ('Q', 'S'): 68, ('Q', 'T'): 42, ('Q', 'W'): 130, ('Q', 'Y'): 99, ('Q', 'V'): 96,
        ('E', 'A'): 107, ('E', 'G'): 98, ('E', 'H'): 40, ('E', 'I'): 134, ('E', 'L'): 138,
        ('E', 'K'): 56, ('E', 'M'): 126, ('E', 'F'): 140, ('E', 'P'): 93, ('E', 'S'): 80,
        ('E', 'T'): 65, ('E', 'W'): 152, ('E', 'Y'): 122, ('E', 'V'): 121, ('G', 'A'): 60,
        ('G', 'H'): 98, ('G', 'I'): 135, ('G', 'L'): 138, ('G', 'K'): 127, ('G', 'M'): 127,
        ('G', 'F'): 153, ('G', 'P'): 42, ('G', 'S'): 56, ('G', 'T'): 59, ('G', 'W'): 184,
        ('G', 'Y'): 147, ('G', 'V'): 109, ('H', 'A'): 86, ('H', 'I'): 94, ('H', 'L'): 99,
        ('H', 'K'): 32, ('H', 'M'): 87, ('H', 'F'): 100, ('H', 'P'): 77, ('H', 'S'): 89,
        ('H', 'T'): 47, ('H', 'W'): 115, ('H', 'Y'): 83, ('H', 'V'): 84, ('I', 'A'): 94,
        ('I', 'L'): 5, ('I', 'K'): 102, ('I', 'M'): 10, ('I', 'F'): 21, ('I', 'P'): 95,
        ('I', 'S'): 142, ('I', 'T'): 89, ('I', 'W'): 61, ('I', 'Y'): 33, ('I', 'V'): 29,
        ('L', 'A'): 96, ('L', 'K'): 107, ('L', 'M'): 15, ('L', 'F'): 22, ('L', 'P'): 98,
        ('L', 'S'): 145, ('L', 'T'): 92, ('L', 'W'): 61, ('L', 'Y'): 36, ('L', 'V'): 32,
        ('K', 'A'): 106, ('K', 'M'): 95, ('K', 'F'): 102, ('K', 'P'): 103, ('K', 'S'): 121,
        ('K', 'T'): 78, ('K', 'W'): 110, ('K', 'Y'): 85, ('K', 'V'): 97, ('M', 'A'): 84,
        ('M', 'F'): 28, ('M', 'P'): 87, ('M', 'S'): 135, ('M', 'T'): 81, ('M', 'W'): 67,
        ('M', 'Y'): 36, ('M', 'V'): 21, ('F', 'A'): 113, ('F', 'P'): 114, ('F', 'S'): 155,
        ('F', 'T'): 103, ('F', 'W'): 40, ('F', 'Y'): 22, ('F', 'V'): 50, ('P', 'A'): 27,
        ('P', 'S'): 74, ('P', 'T'): 38, ('P', 'W'): 147, ('P', 'Y'): 110, ('P', 'V'): 68,
        ('S', 'A'): 99, ('S', 'T'): 58, ('S', 'W'): 177, ('S', 'Y'): 144, ('S', 'V'): 124,
        ('T', 'A'): 58, ('T', 'W'): 128, ('T', 'Y'): 92, ('T', 'V'): 69, ('W', 'A'): 148,
        ('W', 'Y'): 37, ('W', 'V'): 88, ('Y', 'A'): 112, ('Y', 'V'): 55, ('V', 'A'): 64
    }
    for (aa1, aa2), distance in list(grantham.items()):
        grantham[(aa2, aa1)] = distance
        grantham[(aa1, aa1)] = 0
        grantham[(aa2, aa2)] = 0
    return grantham

GRANTHAM_MATRIX = grantham_matrix()

#grantham distance between two amino acid
def grantham_distance(aa1, aa2):
    return GRANTHAM_MATRIX.get((aa1.upper(), aa2.upper()), 0)


# read clustal alignment file
def read_clustal(file_name):
    all_alignments = []
    sequences = {}
    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('CLUSTAL') or not line:
                if sequences:
                    all_alignments.append(sequences)
                    sequences = {}
                continue
            parts = line.split()
            if len(parts) >= 2:
                seq_id, seq = parts[0], ''.join(parts[1:])
                if seq_id not in sequences:
                    sequences[seq_id] = ''
                sequences[seq_id] += seq
    if sequences:
        all_alignments.append(sequences)
    return all_alignments

# check alignment for identical, conserved and variable 
def analyze_alignment(alignment):
    identical = 0
    conserved = 0
    variable = 0
    total = 0
    min_length = min(len(seq) for seq in alignment.values())
    for i in range(min_length):
        column = [seq[i] for seq in alignment.values() if i < len(seq)]
        if '-' in column or len(column) != len(alignment):
            continue
        total += 1
        if len(set(column)) == 1:
            identical += 1
        elif all(grantham_distance(a, b) < 100 for a in column for b in column):
            conserved += 1
        else:
            variable += 1
    if total == 0:
        return 0, 0, 0
    return identical / total, conserved / total, variable / total


#comparison between original and custom alignment 
#plot a bar chart comparing original and custom alignment proportions
def comparison(original_props, custom_props):
    categories = ['Identical', 'Conserved', 'Variable']
    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(categories))
    width = 0.35
    #plot bars for original and custom data with slight offsets
    ax.bar([i - width/2 for i in x], original_props, width, label='Original', alpha=0.7)
    ax.bar([i + width/2 for i in x], custom_props, width, label='Custom', alpha=0.7)
    ax.set_ylabel('Proportion')
    ax.set_title('Alignment Comparison')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend()
    plt.tight_layout()
    plt.show()


# Main
def main():
    original_file = input("enter the name of the original alignment file: ")
    custom_file = input("enter the name of the custom alignment file: ")

    original_alignments = read_clustal(original_file)
    custom_alignments = read_clustal(custom_file)

    original_totals = [0, 0, 0]
    custom_totals = [0, 0, 0]

    print("\nProcessing Original Alignments:")
    for i, alignment in enumerate(original_alignments):
        print(f"\nAlignment {i + 1}:")
        props = analyze_alignment(alignment)
        original_totals = [x + y for x, y in zip(original_totals, props)]

    print("\nProcessing Custom Alignments:")
    for i, alignment in enumerate(custom_alignments):
        print(f"\nAlignment {i + 1}:")
        props = analyze_alignment(alignment)
        custom_totals = [x + y for x, y in zip(custom_totals, props)]

    num_original = len(original_alignments)
    num_custom = len(custom_alignments)
    original_props = [x / num_original for x in original_totals]
    custom_props = [x / num_custom for x in custom_totals]

    print("\nOverall Results:")
    print("\nOriginal Alignments:")
    print(f"Identical: {original_props[0]:.2%}")
    print(f"Conserved: {original_props[1]:.2%}")
    print(f"Variable: {original_props[2]:.2%}")

    print("\nCustom Alignments:")
    print(f"Identical: {custom_props[0]:.2%}")
    print(f"Conserved: {custom_props[1]:.2%}")
    print(f"Variable: {custom_props[2]:.2%}")

    comparison(original_props, custom_props)

if __name__ == "__main__":
    main()