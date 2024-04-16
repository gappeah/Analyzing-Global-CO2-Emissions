# DNA/RNA Sequence Converter and Secondary Structure Predictor

This Python project is designed to convert DNA or RNA sequences into amino acids and predict their secondary structure based on the amino acid sequence.

## Functionality

The program takes a DNA or RNA sequence as input and performs the following steps:

1. **Input**: The user is prompted to enter a DNA or RNA sequence.
2. **Transcription (DNA to RNA)**: If the input sequence is DNA, it is transcribed into RNA by replacing thymine (T) with uracil (U).
3. **Translation (RNA to Amino Acids)**: The transcribed RNA sequence is then translated into a sequence of amino acids using a predefined genetic code dictionary.
4. **Secondary Structure Prediction**: Based on the sequence of amino acids, the program predicts the secondary structure (helix, strand, or coil) using the GOR algorithm.
5. **Output**: The program displays the amino acid sequence and the predicted secondary structure.

## Genetic Code Dictionary

The genetic code dictionary maps RNA codons to amino acids according to the standard genetic code. It covers all possible combinations of RNA codons and their corresponding amino acids, including start and stop codons.

## Secondary Structure Prediction

The program predicts the secondary structure of the protein using the GOR algorithm, which assigns probabilities to each amino acid for being part of a helix (H), strand (E), or coil (C). These probabilities are pre-defined based on empirical data.

## Usage

To use the program:

1. Run the script.
2. Enter a DNA or RNA sequence when prompted.
3. The program will display the corresponding amino acid sequence and predicted secondary structure.

## Example Usage

```python
Enter a DNA or RNA sequence: AUGGCGAAGCUAA
Amino Acid Sequence: MGR*
Predicted Secondary Structure: HHHHHHHHCCC
```

## Note

- Invalid inputs, such as sequences containing characters other than A, C, G, T, or U, will result in an error message prompting the user to enter a valid DNA or RNA sequence.