import os
import sys 


def define_arguments(): # Fetch and format arguments

    if len(sys.argv) != 4: # Check that 4 arguments are given
        print("This script requires 4 inputs: path/to/script, path/to/input, sample_id, path/to/output")
        sys.exit(1) # Exit process and raise a SystemExit exception

    input_file_path = sys.argv[1] # Assign the 2nd argument as the input file path
    sample_id = sys.argv[2] # Assign the 3rd argument as the sample id
    output_dir= sys.argv[3] # Assign the 4th argument as the output directory path

    if not output_dir.endswith('/'): # Ensure the output directory path ends with a slash
        output_dir += '/'

    os.makedirs(output_dir, exist_ok=True) # Create the output directory if it does not already exist

    output_file_name=f'{output_dir}trimmed_{sample_id}.bed' # Format the emitted file path

    return input_file_path, output_file_name


def trim_reads(input_file_path, output_file_name): # Trim reads function

    with open(input_file_path, 'r') as infile, open(output_file_name, 'w') as outfile: # Open the input file (read mode) and output file (write mode)
        for line in infile:
            fields = line.strip().split('\t') # Remove whitespace and split based on tab delimiter
            chromosome_id, start, end, read_id, score, direction = fields[:6] # Define data fields

            start, end = int(start), int(end) # Explicity state start and end as integers

            # Perform trimming
            if direction == '+':
                end = start + 1 # If the strand direction is positive trim the end nucleotide 
            elif direction == '-':
                start = end - 1 # If the strand direction is positive trim the start nucleotide

            # Write the trimmed data to the output file
            outfile.write(f"{chromosome_id}\t{start}\t{end}\t{read_id}\t{score}\t{direction}\n")


def main(): # Script workflow

    input_file_path, output_file_path = define_arguments() # Call define_arguments() function

    trim_reads(input_file_path, output_file_path) # Call trim_reads() function and pass it the args from define_arguments()


if __name__ == '__main__': # Safety check to only execute the script if the file is run directly

    main() # Call main() function to execute script workflow