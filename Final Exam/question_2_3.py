"""
Question 2. What is the length of the longest sequence in the file?

Question 3. What is the length of the shortest sequence in the file?
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

# Sort sequences based on length (longest to shortest)
sorted_fasta = dict(sorted(fasta_dict.items(), key=lambda item: len(item[1]), reverse=True))

# Find the longest sequence
longest_seq_name, longest_seq = max(fasta_dict.items(), key=lambda item: len(item[1]))

# Find the shortest sequence
shortest_seq_name, shortest_seq = min(fasta_dict.items(), key=lambda item: len(item[1]))

# Print results
print(f"Longest sequence: {longest_seq_name}, Length: ({len(longest_seq)} bp)")
print(f"Shortest sequence: {shortest_seq_name}, Length: ({len(shortest_seq)} bp)")
