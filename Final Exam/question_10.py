"""
Question 10. Which one of the following repeats of length 7 has a maximum number of occurrences?

a. TGCGCGC
b. GCGCGCA
c. CGCGCCG
d. CATCGCC
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

# Find the most frequent repeat of length 7 or sevenmer in all sequences

from collections import defaultdict

# Dictionary to store sevenmer (7-mer) counts
sevenmer_counts = defaultdict(int)

# Extract all 7-mers from every sequence in fasta_dict
for seq_name, sequence in fasta_dict.items():
    for i in range(len(sequence) - 6):  # Ensure a valid 7-mer
        sevenmer = sequence[i:i+7]  # Extract 7-mer
        sevenmer_counts[sevenmer] += 1  # Count occurrences

# Find the sevenmer with the highest frequency
most_frequent_sevenmer = max(sevenmer_counts, key=sevenmer_counts.get)

print("Most frequent sevenmer repeat: {most_frequent_sevenmer}")
