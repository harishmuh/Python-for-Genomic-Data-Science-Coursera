"""
Question 4. What is the length of the longest ORF appearing in reading frame 2 of any of the sequences?
"""

# Function for opening a fasta file 
def open_fasta_file(filename):
    """
    Opens a FASTA file and reads its contents.
    
    Parameters:
        filename (str): Path to the FASTA file.
    
    Returns:
        list: A list of lines from the file if successfully read.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        IOError: If an error occurs while reading the file.
    """
    try:
        with open(filename, "r") as file:
            return file.readlines()  # Read all lines and return as a list  
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None
    except IOError:
        print(f"Error: Could not read the file '{filename}'.")
        return None

# Function to parse the fasta file and store sequences in a dictionary 
def parse_fasta(file):
    """
    Parses a FASTA file and stores sequences in a dictionary.

    Args:
        file (file object): Opened FASTA file.

    Returns:
        dict: A dictionary with sequence names as keys and sequences as values.
    """
    sequences = {}  # Dictionary to store sequences
    name = None  # Initialize sequence name

    for line in file:
        line = line.rstrip()  # Remove trailing newline characters

        if line.startswith('>'):  # Identify sequence header
            words = line.split()  # Split header line into words
            name = words[0][1:]  # Extract sequence name (without '>')
            sequences[name] = ''  # Initialize sequence in dictionary
        else:
            if name:  # Ensure there's an active sequence name
                sequences[name] += line  # Append sequence data

    return sequences

# Parsing sequences into dictionary
filename = "dna2.fasta"
file = open_fasta_file(filename)
fasta_dict = parse_fasta(file)

# Creating function for finding ORFs in a reading frame
def find_orfs_in_frame(sequence, frame=1): # 1 as default frame 
    """
    Finds Open Reading Frames (ORFs) in a given reading frame.

    Parameters:
        sequence (str): The DNA sequence (uppercase).
        frame (int): The reading frame (1, 2, or 3).

    Returns:
        list: A list of ORFs found in the given reading frame.
    """
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []

    seq = sequence[frame - 1:]  # Adjust reading frame
    seq_length = len(seq)

    i = 0
    while i < seq_length - 2:
        codon = seq[i:i+3]

        if codon == start_codon:
            start_index = i
            for j in range(i, seq_length - 2, 3):
                stop_codon = seq[j:j+3]
                if stop_codon in stop_codons:
                    orfs.append(seq[start_index:j+3])  # Include stop codon
                    i = j  # Move to next possible ORF
                    break
        i += 3  # Move to next codon

    return orfs

# Function to find ORFs in all sequences 
def find_orfs_in_fasta(fasta_dict):
    """
    Finds ORFs in each reading frame (1, 2, 3) for all sequences in a FASTA dictionary.

    Parameters:
        fasta_dict (dict): Dictionary with sequence names as keys and DNA sequences as values.

    Returns:
        dict: Dictionary with sequence names as keys and a sub-dictionary of ORFs in each frame.
    """
    orf_results = {}

    for name, sequence in fasta_dict.items():
        sequence = sequence.upper()  # Ensure uppercase for consistency
        orf_results[name] = {
            "Frame 1": find_orfs_in_frame(sequence, frame=1),
            "Frame 2": find_orfs_in_frame(sequence, frame=2),
            "Frame 3": find_orfs_in_frame(sequence, frame=3)
        }

    return orf_results

# Creating dictionary that contains all ORFs in all frame from all sequences
orf_dict = find_orfs_in_fasta(fasta_dict)

# Function to get the longest ORF in specific frame

def get_longest_orf_by_frame(orf_dict, frame):
    """
    Finds the longest ORF from a specified reading frame across all sequences.

    Parameters:
        orf_dict (dict): Dictionary containing ORFs for each sequence.
        frame (str): Reading frame to search in ("Frame 1", "Frame 2", or "Frame 3").

    Returns:
        tuple: (Longest ORF sequence, Length of the longest ORF, Sequence name)
    """
    longest_orf = ""
    longest_length = 0
    longest_seq_name = ""

    for name, frames in orf_dict.items():
        if frame in frames:  # Ensure the requested frame exists in the dictionary
            for orf in frames[frame]:  # Iterate through all ORFs in the specified frame
                if len(orf) > longest_length:
                    longest_orf = orf
                    longest_length = len(orf)
                    longest_seq_name = name

    return longest_seq_name, longest_length, longest_orf

# The longest ORF in Frame 2
frame_to_check_fr2 = "Frame 2"
seq_name_fr2, longest_orf_length_fr2, longest_orf_seq_fr2 = get_longest_orf_by_frame(orf_dict, frame_to_check_fr2)

# Print result
print(f"Longest ORF in {frame_to_check_fr2}:")
print(f"Sequence Name: {seq_name_fr2}")
print(f"Length: {longest_orf_length_fr2} bp")
