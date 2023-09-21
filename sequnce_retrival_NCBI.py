from Bio import Entrez

def fetch_sequence(accession_number):
    Entrez.email = "nishantthalwal@gmail.com"  # Enter your email address here
    try:
        # Fetch the sequence using the accession number
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="fasta", retmode="text")
        sequence = handle.read()
        handle.close()
        return sequence
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    accession_number = input("Input The Accession Number")  # Replace with the accession number you want to retrieve
    sequence = fetch_sequence(accession_number)
    
    if sequence:
        with open(f"{accession_number}.fasta", "w") as output_file:
            output_file.write(sequence)
        print(f"Sequence saved to {accession_number}.fasta")
