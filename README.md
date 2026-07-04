# Open Reading Frames Assignment
Find the longest open reading frame from a group of fasta DNA sequences

## Usage
1. Make sure you have Python installed on your system.
2. Run the script by providing the input file as a command-line argument. Optionally, include `--print_length` to print the length:

    ```bash
    python orf_finder.py input_file.fasta --print_length
    ```

## Input File Format
The input file will be in the fasta format, containing DNA sequences.

## Output
The script will print the longest DNA ORF found in the input sequences. If `--print_length` is specified, it will also print the length.

## Example
```bash
# Without printing length
python orf_finder.py example.txt

# With printing length
python orf_finder.py example.txt --print_length
```
