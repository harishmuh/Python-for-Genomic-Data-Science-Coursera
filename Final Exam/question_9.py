"""
Question 9. Find all repeats of length 12 in the input file. Let's use Max to specify the number of copies of the most frequent repeat of length 12. 
How many different 12-base sequences occur Max times?
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

# Find all repeat of length 12 or twelevemer

# Importing libary 'defaultdict'
from collections import defaultdict

# Dictionary to store 12-mer (12-base) counts
twelvemer_counts = defaultdict(int)

# Extract all 12-mers from every sequence in fasta_dict
for seq_name, sequence in fasta_dict.items():
    for i in range(len(sequence) - 11):  # Ensure valid 12-mer
        twelvemer = sequence[i:i+12]  # Extract 12-mer
        twelvemer_counts[twelvemer] += 1  # Count occurrences

# Find the maximum frequency (Max)
max_frequency = max(twelvemer_counts.values())

# Count how many 12-mers occur exactly 'Max' times
num_max_occurrences = sum(1 for count in twelvemer_counts.values() if count == max_frequency)

# Print results
print(f"The most frequent 12-mer (Max frequency): {max_frequency}")
print(f"Number of different 12-base sequences that occur {max_frequency} times: {num_max_occurrences}")
