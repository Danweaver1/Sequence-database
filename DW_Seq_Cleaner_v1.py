import sys
from Bio import SeqIO

### Modified version of "Seq_cleaner.py" from Biopython (genivaldo.gueiros@gmail.com)
# Modifications:
#	Added a max sequence length parameter - only sequences equal to or less than this value are kept
#	Output a file of sequences which fail QC.
#	Output a file of sequences which weren't added due to being duplicates.
# To run: 
#	python DW_Seq_Cleaner.py <fasta file> <max_length> <min_length> <min percent Ns>

def sequence_cleaner(fasta_file, max_length=0, min_length=0, por_nucl=0):
    # Create our hash table to add the sequences
    sequences={}
    QCfail={}
    dupl={}

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) <= max_length and len(sequence) >= min_length and
            (float(sequence.count("A")) + float(sequence.count("C")) + float(sequence.count("G")) + float(sequence.count("T")))/float(len(sequence)))*100 >= por_nucl:
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its header are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.description
       # If it is already in the hash table, we're just going to add to dupl
            else:
                dupl[sequence] = seq_record.description
       #any that fail parameters go into QC fail dictionary
        else:
            QCfail[sequence] = seq_record.description     


    # Write the clean sequences

    # Create a file in the same directory where you ran this script
    with open("clear_" + fasta_file, "w+") as output_file:
        # Just read the sequence dictionary and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")
# Write the sequences whcih failed QC
    with open("QCfail_" + fasta_file, "w+") as output_file:
        # Just read the QC fail dictionary and write on the file as a fasta format
        for sequence in QCfail:
            output_file.write(">" + QCfail[sequence] + "\n" + sequence + "\n")
# Write the duplicated sequences
    with open("duplicates_" + fasta_file, "w+") as output_file:
        # Just read the duplicates dictionary and write on the file as a fasta format
        for sequence in dupl:
            output_file.write(">" + dupl[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + fasta_file)


userParameters = sys.argv[1:]
# check enough command line parameters have been given
try:
    if len(userParameters) == 1:
        sequence_cleaner(userParameters[0])
    elif len(userParameters) == 2:
        sequence_cleaner(userParameters[0], float(userParameters[1]))
    elif len(userParameters) == 3:
        sequence_cleaner(userParameters[0], float(userParameters[1]),
                         float(userParameters[2]))
    elif len(userParameters) == 4:
        sequence_cleaner(userParameters[0], float(userParameters[1]),
                         float(userParameters[2]),float(userParameters[3]))
    else:
        print("There is a problem!")
except:
    print("There is a problem!")

