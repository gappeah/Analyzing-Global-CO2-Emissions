# Define the genetic code dictionary
genetic_code = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
    'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W'
}

from Bio.PDB import PDBParser, PPBuilder
import nglview as nv
from io import StringIO

def visualize_secondary_structure(amino_acid_sequence):
    # Create a string representation of the amino acid sequence
    sequence_str = ''.join(amino_acid_sequence)

    # Use Biopython to create a protein structure
    structure = PPBuilder().build_peptide(sequence_str)

    # Write the structure to a StringIO object
    structure_io = StringIO()
    io = PDBIO()
    io.set_structure(structure)
    io.save(structure_io)

    # Create an NGLView instance and load the structure
    view = nv.show_structure_file(structure_io.getvalue())

    # Add secondary structure representation
    view.add_representation("cartoon", selection="all", color="residueindex")

    # Display the structure
    view
    
def get_input():
    sequence = input("Enter a DNA or RNA sequence: ")
    return sequence.upper()

def transcribe_dna_to_rna(dna):
    return dna.replace('T', 'U')

def translate_rna_to_amino_acids(rna):
    amino_acids = []
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if codon in genetic_code:
            amino_acid = genetic_code[codon]
            amino_acids.append(amino_acid)
        else:
            break
    return amino_acids

def predict_secondary_structure(amino_acid_sequence):
    # GOR probability parameters
    gor_helix = [0.888, 0.544, 1.032, -0.032, 1.005, 0.285, 0.954, 0.315, -0.384, -0.284, -0.695, -0.072, -0.647, -0.192, -1.612, -0.232, -0.368, -1.036, -0.293, -0.057]
    gor_strand = [0.836, -0.920, -0.179, 1.603, -0.590, -1.892, -0.109, -1.649, -0.329, -1.263, -1.053, -1.370, -1.066, -0.621, -1.022, -0.506, -0.355, -0.766, -1.075, -0.697]
    gor_coil = [-0.728, -0.359, -0.295, -0.543, -0.292, 0.307, -0.365, 0.417, 0.603, 0.570, 1.212, 1.118, 1.185, 0.441, 1.924, 0.045, 0.194, 1.170, 0.933, 0.238]

    # Initialize lists to store predicted structures
    helix_prob = []
    strand_prob = []
    coil_prob = []

    # Iterate through the amino acid sequence
    for amino_acid in amino_acid_sequence:
        if amino_acid == '*':  # Handle stop codons separately
            helix_prob.append(0)
            strand_prob.append(0)
            coil_prob.append(0)
        else:
            # Get the index of the amino acid in the probability lists
            index = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'].index(amino_acid)

            # Append the probabilities to the respective lists
            helix_prob.append(gor_helix[index])
            strand_prob.append(gor_strand[index])
            coil_prob.append(gor_coil[index])

    # Initialize lists to store predicted structures
    predicted_structure = []

    # Iterate through the probabilities and predict the structure
    for i in range(len(amino_acid_sequence)):
        max_prob = max(helix_prob[i], strand_prob[i], coil_prob[i])
        if max_prob == helix_prob[i]:
            predicted_structure.append('H')
        elif max_prob == strand_prob[i]:
            predicted_structure.append('E')
        else:
            predicted_structure.append('C')

    return ''.join(predicted_structure)

def main():
    sequence = get_input()
    if set(sequence) <= set('ACGT'):
        rna = transcribe_dna_to_rna(sequence)
        amino_acids = translate_rna_to_amino_acids(rna)
        secondary_structure = predict_secondary_structure(amino_acids)
        print("Amino Acid Sequence:", ''.join(amino_acids))
        print("Predicted Secondary Structure:", secondary_structure)
    elif set(sequence) <= set('ACGU'):
        amino_acids = translate_rna_to_amino_acids(sequence)
        secondary_structure = predict_secondary_structure(amino_acids)
        print("Amino Acid Sequence:", ''.join(amino_acids))
        print("Predicted Secondary Structure:", secondary_structure)
    else:
        print("Invalid input. Please enter a valid DNA or RNA sequence.")

if __name__ == "__main__":
    main()