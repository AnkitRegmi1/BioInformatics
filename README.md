# Bioinformatics Final Assignment — Custom Substitution Matrix, Progressive MSA, and Alignment Evaluation

This repository contains the code and data for a bioinformatics final assignment. The project:

- Builds a custom amino‑acid substitution matrix (log‑odds) from curated multiple alignments of globin families.
- Performs a progressive multiple sequence alignment (MSA) of ABHD11 protein homologs using the custom matrix.
- Compares “original” vs “customized” alignments and evaluates column quality (identical, conserved, variable) using Grantham distances, with visualization.

## Repository structure

- `grouped_alignments.txt`: Curated grouped alignments (e.g., Cytoglobin, Hemoglobin α/β, Myoglobin, Neuroglobin) used to derive pairwise statistics and the custom matrix.
- `test.py`: Processes `grouped_alignments.txt` to compute:
  - Pairwise occurrence counts
  - Observed and expected pair probabilities
  - Log‑odds substitution scores (saved to CSV)
  - Per‑group and total amino‑acid frequencies
- Outputs produced by `test.py`:
  - `combined_occurrences_matrix.csv`
  - `combined_observed_probability_matrix.csv`
  - `combined_expected_probability_matrix.csv`
  - `combined_log_odds_matrix.csv` (used as substitution matrix downstream)
  - `total_amino_acid_frequencies.csv` and per‑group frequency CSVs
- `sequence.txt`: FASTA‑like file with ABHD11 protein sequences (multiple species) used for progressive MSA.
- `question3.py`: Implements progressive MSA using:
  - The custom substitution matrix in `combined_log_odds_matrix.csv`
  - A fixed gap penalty
  - Outputs a formatted alignment to stdout
- `originalalignment.txt`, `customizedalignment.txt`: Blocked alignments for comparison (original vs customized) used by the evaluator.
- `milestone3.py`: Evaluates alignments by column using Grantham distances:
  - Identical (all residues equal)
  - Conserved (all pairwise Grantham < 100)
  - Variable (otherwise)
  - Prints proportions and shows a bar chart comparison.
- Additional data files:
  - `sequences.fasta`: Example ABHD2 homologs (not directly used by `question3.py`).
  - `ABHD*.txt`: Additional sequence sets used during exploration.

## Requirements

- Python 3.8+ (tested with standard library modules only, except plotting)
- matplotlib (for `milestone3.py` visualization)

Install matplotlib:

```bash
pip install matplotlib
```

## Reproducible workflow

1) Generate the custom substitution matrix from grouped alignments

```bash
python test.py
```

This reads `grouped_alignments.txt` and writes the CSV outputs listed above, most importantly `combined_log_odds_matrix.csv`.

2) Run progressive MSA on ABHD11 sequences using the custom matrix

```bash
python question3.py
```

What it does:
- Reads `combined_log_odds_matrix.csv` and sequences from `sequence.txt`
- Performs progressive alignment (global alignment against a growing reference)
- Prints a formatted multiple alignment to the console

Notes:
- The gap penalty is set in code (currently `-50`), chosen to be larger in magnitude than any single entry in the log‑odds matrix.
- `sequence.txt` expects FASTA headers starting with `>` and uses the first whitespace‑separated token as the sequence ID.

3) Evaluate “original” vs “customized” alignments with Grantham‑based categories

```bash
python milestone3.py
```

When prompted, provide the two alignment files, for example:
- Original: `originalalignment.txt`
- Customized: `customizedalignment.txt`

What it does:
- Parses block‑style alignments (sequence ID followed by aligned residues; blocks separated by blank lines)
- Classifies each ungapped column as identical, conserved (by Grantham < 100), or variable
- Prints overall proportions and displays a bar chart comparison

## Data formats

- Grouped alignments (`grouped_alignments.txt`):
  - Headings end with `:` (e.g., `Cytoglobin A:`) followed by sequence lines; blank lines separate blocks.
- FASTA‑like sequences (`sequence.txt`):
  - `>` header lines, sequence ID is the first token on the header line; sequence lines may wrap.
- Alignment comparison files (`originalalignment*.txt`, `customizedalignment*.txt`):
  - Each block contains lines of `SEQ_ID<spaces>ALIGNED_RESIDUES`, separated by blank lines.

## Key implementation details

- Custom matrix construction (`test.py`):
  - Counts pair co‑occurrences per column across normalized sequence blocks (ignores gaps, applies +0.5 to upper triangle)
  - Computes observed probabilities q(ij), expected probabilities e(ij) (from diagonal marginals), and log‑odds scores s(ij) = 2·round(log2(q/e))
  - Saves in upper‑triangular CSV format with amino‑acid headers
- Progressive MSA (`question3.py`):
  - Needleman–Wunsch style global alignment with the custom substitution scores and fixed gap penalty
  - Progressively aligns each sequence to the evolving reference and normalizes lengths
- Alignment evaluation (`milestone3.py`):
  - Uses an explicit Grantham distance map; columns with all pairwise distances < 100 are “conserved”

## Example commands (Windows PowerShell or any shell)

```bash
# 1) Build the matrix and derived CSVs
python test.py

# 2) Progressive MSA of ABHD11 sequences
python question3.py

# 3) Evaluate and visualize original vs customized alignments
python milestone3.py
# When prompted, enter:
#   enter the name of the original alignment file: originalalignment.txt
#   enter the name of the custom alignment file: customizedalignment.txt
```

## Results and outputs

- CSV matrices in the repository root after step (1)
- A printed multiple alignment after step (2)
- Console summary and a matplotlib bar chart after step (3)

## Acknowledgments

- Grantham, R. (1974). Amino acid difference formula to help explain protein evolution.
- Course staff and assignment specification for guidance and datasets.

## License

No license specified. If you intend to reuse this work, please add an appropriate license.

