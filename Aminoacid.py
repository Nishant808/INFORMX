from Bio import SeqIO
from Bio.Seq import Seq

# Function to read a FASTA file
def read_fasta_file(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

# Function to translate a DNA or RNA sequence to amino acids
def translate_to_amino_acids(dna_sequence, genetic_code="Standard"):
    try:
        sequence = Seq(dna_sequence)
        protein_sequence = sequence.translate(table=genetic_code)
        return protein_sequence
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Example usage
if __name__ == "__main__":
    # Specify the path to your FASTA file containing DNA or RNA sequence
    fasta_file_path = "your_sequence.fasta"

    # Read the FASTA file
    sequences = read_fasta_file(fasta_file_path)

    # Assuming you have a single sequence in the file, you can access it like this:
    dna_sequence = sequences[0].seq

    # Translate the DNA or RNA sequence to amino acids (protein sequence)
    protein_sequence = translate_to_amino_acids(dna_sequence)

    if protein_sequence:
        print("Amino Acid Sequence:")
        print(protein_sequence)
    else:
        print("Translation failed.")
