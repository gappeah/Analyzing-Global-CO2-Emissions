# DNA/RNA Sequence Converter and Secondary Structure Predictor

This Python script allows users to convert DNA or RNA sequences, predict the secondary structure of proteins, and output the results to a CSV file.

## Dependencies
- Python 3.x
- `csv` module
- `input()` function for user input

## Installation
No installation is required. Simply download the script and run it using Python.

## Usage
1. Run the script.
2. Enter a DNA or RNA sequence when prompted.
3. The script will:
   - Convert the input sequence to uppercase.
   - Transcribe DNA sequences to RNA.
   - Translate RNA sequences to amino acids.
   - Predict the secondary structure of the protein using the GOR method.
   - Output the amino acid sequence and predicted secondary structure to a CSV file named `output.csv`.

## Example
```
Enter a DNA or RNA sequence: AUGCUCAAGU
```
Output:
```
Amino Acid Sequence,Predicted Secondary Structure
MLQQH
CCCEE
```

## Script Explanation
- `genetic_code`: Dictionary containing the genetic code.
- `get_input()`: Function to prompt the user for a DNA or RNA sequence input.
- `transcribe_dna_to_rna(dna)`: Function to transcribe DNA sequences to RNA.
- `translate_rna_to_amino_acids(rna)`: Function to translate RNA sequences to amino acids.
- `predict_secondary_structure(amino_acid_sequence)`: Function to predict the secondary structure of proteins using the GOR method.
- `main()`: Main function to execute the sequence conversion and secondary structure prediction, and output the results to a CSV file.

## Output Format
The output CSV file contains two columns:
1. Amino Acid Sequence: The translated sequence of amino acids.
2. Predicted Secondary Structure: The predicted secondary structure corresponding to each amino acid.