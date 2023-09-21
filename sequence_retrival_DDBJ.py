import requests

def retrieve_sequence_from_ddbj(accession_number):
    # Define the DDBJ API URL
    ddbj_api_url = f"https://getentry.ddbj.nig.ac.jp/rest/entry/nucleotide/{accession_number}.fasta"

    try:
        # Send a GET request to the DDBJ API
        response = requests.get(ddbj_api_url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # The sequence data is in the response content
            sequence_data = response.text
            return sequence_data
        else:
            print(f"Failed to retrieve sequence. Status code: {response.status_code}")
            return None
    except requests.exceptions.RequestException as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    accession_number = input("Input The Accession Number")  # Replace with the accession number you want to retrieve
    sequence = retrieve_sequence_from_ddbj(accession_number)
    if sequence:
        print(f"Sequence for accession {accession_number}:\n{sequence}")
