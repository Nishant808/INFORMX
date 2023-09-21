from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Function to read a FASTA file
def read_fasta_file(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences

# Function to perform BLAST alignment
def perform_blast_alignment(query_sequence, database_name="nr", evalue=0.001):
    try:
        # Perform BLAST search
        result_handle = NCBIWWW.qblast("blastn", database_name, query_sequence, expect=evalue)

        # Parse BLAST result
        blast_records = NCBIXML.parse(result_handle)
        return blast_records
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

# Example usage
if __name__ == "__main__":
    # Specify the path to your FASTA file
    fasta_file_path = "your_sequence.fasta"

    # Read the FASTA file
    sequences = read_fasta_file(fasta_file_path)

    # Assuming you have a single query sequence, you can access it like this:
    query_sequence = sequences[0]

    # Perform BLAST alignment
    blast_records = perform_blast_alignment(query_sequence.seq)

    if blast_records:
        # Print alignment results
        for record in blast_records:
            print(f"Alignment for {record.query_id}:")
            for alignment in record.alignments:
                print(f"Subject: {alignment.title}")
                for hsp in alignment.hsps:
                    print(f"Length: {hsp.align_length}")
                    print(f"Score: {hsp.score}")
                    print(f"E-value: {hsp.expect}")
                    print(f"Alignment: {hsp.query}\n{hsp.match}\n{hsp.sbjct}\n")
    else:
        print("No BLAST records found.")
