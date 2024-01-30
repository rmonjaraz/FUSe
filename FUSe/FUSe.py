#!/usr/bin/env python
# coding: utf-8
__version__="v2.2"

"""
@author: Rodrigo Monjaraz-Ruedas (monroderik@gmail.com)

This script performs a series of operations to align, trim and filter alignments.

Alignment is performed using MAFFT as inplemented in  phyluce_align_seqcap_align.

Trimming of alignments is performed with trimAL or Gblocks.

Filtering of sequences and alignments is as follows:

Sequences:
    Remove Short - Removes short sequences based on a percentage of gaps acording with entire lenght of the alignment
    Remove Divergent - Removes divergent sequences based on a pairwise identity comparisson

Alignments:
    No. taxa - Remove alignments that doesn't meet the minimum number of taxa desired
    Alignment lenght - Remove alignments that doesn't meet the minimum lenght desired (in bp)
    Completeness - Generates completenes matricess based on a percentage of taxa (Similar to min_taxa but based on total number of taxa)


This workflow uses:
Biopython
phyluce: "phyluce_align_seqcap_align" 
numpy
Python 3.6 or higher

make sure you have phyluce installed and cite it accordingly.

Check README.md for detailed instructions on installation, useage and examples.

"""

import re
import os
import sys
import math
import glob
import shutil
import logging
import argparse
import tempfile
import subprocess
import numpy as np
from Bio.Seq import Seq
from datetime import datetime
from Bio import AlignIO, SeqIO
from collections import Counter
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def get_args():
    parser = argparse.ArgumentParser('FUSe',
        description="""Align, Trimm and Filter UCE Sequences and Alignments.""")

    parser_group = parser.add_argument_group("Required arguments")

    parser_group.add_argument(
        "-i", "--input",
        dest="monolithic_file",
        required=True,
        help="""The input monolithic fasta file."""
        )
    parser_group.add_argument(
        "-t","--taxa",
        required=True,
        type=int,
        help="The total number of taxa in all alignments.",
    )
    parser.add_argument(
        "-p", "--prefix",
        default="OUTPUT",
        help="""The name to be used for all output folders (default: OUTPUT).""",
    )
    parser.add_argument(
        "-o", "--out-format",
        choices=["fasta", "nexus", "phylip", "clustal", "emboss"],
        default="fasta",
        help="""The ouput alignment format (default: fasta)."""
        )
    parser.add_argument(
        "-c", "--cores",
        type=int,
        default=1,
        help="""The number of PHYSICAL CPUs (default: 1).""",
    )
    parser.add_argument(
        '--trimAL',
        action='store_true',
        help="""Wether to trim alignments using trimal.""",
    )
    parser.add_argument(
        "-a", "--t-method",
        choices=["automated1", "gappyout", "strict", "strictplus", "nogaps"],
        default="automated1",
        help="""trimAl automated method to use (default: automated1).""",
    )
    parser.add_argument(
        '--gblocks',
        action='store_true',
        help="""Wether to trim alignments using gblocks.""",
    ) 
    parser.add_argument(
        "--b1",
        type=float,
        default=0.5,
        help="""The GBLOCKS -b1 proportion (default: 0.5)""",
    )
    parser.add_argument(
        "--b2",
        type=float,
        default=0.70,
        help="""The GBLOCKS -b2 proportion (default: 0.70)""",
    )
    parser.add_argument(
        "--b3",
        type=int,
        default=10,
        help="""The GBLOCKS -b3 integer value (default: 10)""",
    )
    parser.add_argument(
        "--b4",
        type=int,
        default=4,
        help="""The GBLOCKS -b4 integer value (default: 4)""",
    )   
    parser.add_argument(
        '--remove-div',
        action='store_true',
        help="""Wether to remove divergent sequences from alignments.""",
    )
    parser.add_argument(
        "-d", "--divergent",
        default=0.7,
        type=float,
        help="""Percentage of pairwise identity in every sequence to retain (default: 0.7).""",
    )
    parser.add_argument(
        '--remove-short',
        action='store_true',
        help="""Wether to remove short sequences from alignments.""",
    )
    parser.add_argument(
        "-s", "--short-cutoff",
        default=0.7,
        type=float,
        help="""Percentage of "-" (gaps) in every sequece to retain (default: 0.7).""",
    )
    parser.add_argument(
        '--filter-alignments',
        action='store_true',
        help="""Wether to filter alignments by No. taxa and length in bp.""",
    )
    parser.add_argument(
        "-l", "--min-length",
        default=50,
        type=int,
        help="""The minimum alignment lenght to retain alignment (default: 50).""",
    )
    parser.add_argument(
        "-m", "--min-taxa",
        default=4,
        type=int,
        help="""The minimum number of taxa in alignments to retain (default: 4).""",
    )
    parser.add_argument(
        '--get-completeness',
        action='store_true',
        help="""Wether to filter alignments by completeness.""",
    )
    parser.add_argument(
        "-e", "--percent",
        default=0.8,
        type=float,
        help="""Completeness matrix percentage to output (default: 0.8).""",
    )
    parser.add_argument(
        '--taxa-count',
        action='store_true',
        help="""Wether to print taxa count in alignments summary.""",
    )
    parser.add_argument(
        '-v', '--version', 
        action='version', 
        version='%(prog)s ' + __version__
    )

    if len(sys.argv) == 1:
        print(__doc__)
        parser.print_help()
        parser.exit(1)

    return parser.parse_args()

# Steup Logger
def setup_logging(args):
    program = "FUSe"
    logger = logging.getLogger(program)

    #Set File and console loggers (to print in both)
    console = logging.StreamHandler()
    logfile = logging.FileHandler(f"{program}.log")
        
    # Set leveles of info for logger
    logger.setLevel(logging.INFO)
    console.setLevel(logging.INFO)
    logfile.setLevel(logging.INFO)
        
    formatter = logging.Formatter(
        "%(levelname)s - %(message)s"
    )
    console.setFormatter(formatter)
    logfile.setFormatter(formatter)
    logger.addHandler(console)
    logger.addHandler(logfile)

    # Text Here will be printed at biggining of the Log file
    text = f" Starting {program} "
    logger.info(text.center(55, "="))

    now = datetime.now()
    dt_now = now.strftime("%d/%m/%Y %H:%M:%S")
    text = f" Date {dt_now} "
    logger.info(text.center(55, "+"))

    # Pass all the arguments used
    for arg, value in sorted(vars(args).items()):
        logger.info(f"Argument --{arg}: {value}")
    
    return logger

# Function to handle directories
def make_dir(log, directory):
    d = os.path.abspath(directory)
    if os.path.exists(d):
        log.critical(f"Output directory {d} exists, "
            f"provide a new prefix or change Output folder name")
        sys.exit()
    os.mkdir(d)
    return d

# Trim alignments with trimAL
def trimal_worker(args):
    in_file, ucename, output_dir, b1, b2, b3, b4, t_method = args
    cmd = (f'trimal -in {in_file} -out {output_dir}/{ucename}.fasta -{t_method} -fasta')
    subprocess.run(cmd, shell=True, check=True)

# Trim alignments with gblocks
def gblocks_worker(args):
    in_file, ucename, output_dir, b1, b2, b3, b4, t_method = args
    # determine the number of sequences in align
    with open(in_file, "r") as input_aln:
        aln = AlignIO.read(input_aln, "fasta")
        taxa = len(aln)
    b1 = int(round(b1 * taxa)) + 1
    b2 = int(round(b2 * taxa))
    if b2 < b1:
        b2 = b1
    
    cmd = (f' Gblocks {in_file} -t=DNA -b1={b1} -b2={b2} '
        f'-b3={b3} -b4={b4} -b5=h -p=n')
    subprocess.run(cmd, shell=True, check=False)
    
    # Rename and move files
    trimmed_aln = f'{in_file}-gb'
    output_file = os.path.join(output_dir, f'{ucename}.fasta')
    os.rename(trimmed_aln, output_file)


def trim_aligns(log, function, in_dir, output_dir, t_method, b1, b2, b3, b4, prefix, cores):
    output_dir_trimmed = make_dir(log, f"{output_dir}/{prefix}-02-trimmed")
    alignments_list = glob.glob(os.path.join(in_dir, "*.fasta"))    
    args_list = [(fasta_file, 
                os.path.splitext(os.path.basename(fasta_file))[0], 
                output_dir_trimmed, 
                b1, b2, b3, b4, 
                t_method) 
                for fasta_file in alignments_list]
    with Pool(cores) as p:
        p.map(function, args_list)

    return output_dir_trimmed

def remove_locus(args):
    file_list, output_dir = args
    for fasta_file in file_list:
        # Parse sequences in each alignment
        alignments = AlignIO.read(fasta_file, 'fasta')
        ucename = os.path.splitext(os.path.basename(fasta_file))[0].strip("_phased")
        new_align = MultipleSeqAlignment([])
        for record in alignments:
            new_seq_name = re.sub("^(_R_)*{}_*".format(ucename), "", record.id)
            record.id = new_seq_name
            record.description = new_seq_name
            new_align.append(record)
    
        # Write the new alignment to a single output file
        output_file = os.path.join(output_dir, f'{ucename}.fasta')
        with open(output_file, "w") as outf:
            AlignIO.write(new_align, outf, "fasta")
    sys.stdout.write(".")
    sys.stdout.flush()

# Multiprocessing rename files
def batch_rename(log, input_dir, output_dir, prefix, cores):
    file_list = glob.glob(os.path.join(input_dir, "*.fasta"))
    output_dir_align = make_dir(log, f"{output_dir}/{prefix}-01-mafft")
    args = [(file_list, output_dir_align) for file in file_list]
    with Pool(cores) as p:
        p.map(remove_locus, args)
    return output_dir_align


# Convert Alignmets
def convert_alignments(log, input_dir, input_format, output_dir, output_format, prefix):
    output_dir_conv = make_dir(log, f"{output_dir}/{prefix}-Final-Alignments-{output_format}")
    for alignment in glob.glob(os.path.join(input_dir, "*.fasta")):
        fname = os.path.splitext(os.path.basename(alignment))[0]
        output_file = f"{output_dir_conv}/{fname}.{output_format}"
        AlignIO.convert(alignment, input_format, output_file, 
            output_format, molecule_type= "DNA")
    return output_dir_conv

# Remove Short Sequences
def fasta_drop(log, alignments_dir, output_dir, drop_cutoff, prefix):
    output_folder = make_dir(log, f"{output_dir}/{prefix}-04-short-removed")

    for input_file in glob.glob(os.path.join(alignments_dir, "*.fasta")):
        ucename = os.path.splitext(os.path.basename(input_file))[0]
        output_filename = f"{output_folder}/{ucename}.fasta"
        
        with open(input_file, "r") as infile, open(output_filename, "w") as droppedfile:
            for seqs in SeqIO.parse(infile, "fasta"):
                name = seqs.id
                seq = seqs.seq
                seqLen = len(seqs)
                gap_count = 0
                for z in range(seqLen):
                    if seq[z] == "-":
                        gap_count += 1
                if (gap_count / float(seqLen)) >= drop_cutoff:
                    # Write Report
                    log.info(f"{name} from {ucename}")
                    #log.warn(f"{name} from {input_file} was removed.")
                else:
                    SeqIO.write(seqs, droppedfile, "fasta")
    return output_folder

# Function to calculate similarity for Divergent Sequences
def calculate_similarity(args):
    seq1, sequences = args
    matches = 0
    total_positions = 0

    for j, base in enumerate(seq1):
        if base != "-":
            other_bases = sequences[:, j]
            other_bases = other_bases[other_bases != "-"]
            counts = np.unique(other_bases, return_counts=True)

            if counts[1].size > 0:
                most_common_count = max(counts[1])
                most_common_bases = counts[0][counts[1] == most_common_count]

                if base in most_common_bases:
                    matches += 1

            total_positions += 1

    similarity = (matches / total_positions) if total_positions > 0 else 0
    return similarity

# Remove Gaps Only
def remove_gap_columns(args):
    output_dir, prefix, fasta_file = args
    #output_dir = make_dir(log, f"{prefix}-alignments-divergent-removed")

    ucename = os.path.splitext(os.path.basename(fasta_file))[0]

    # Parse sequences in each alignment
    alignments = AlignIO.read(fasta_file, "fasta")

    # Identify columns with only gaps
    gap_columns = set()
    for column in range(alignments.get_alignment_length()):
        if all(seq_record.seq[column] == "-" for seq_record in alignments):
            gap_columns.add(column)

    # Create a new alignment without gap columns
    new_alignments = []
    for seq_record in alignments:
        new_seq = Seq("".join(seq_record.seq[i] for i in range(alignments.get_alignment_length()) if i not in gap_columns))
        new_record = SeqRecord(new_seq, id=seq_record.id, description=seq_record.description)
        new_alignments.append(new_record)

    # Write the new alignment to a file
    output_file = os.path.join(output_dir, f"{ucename}.fasta")
    AlignIO.write(MultipleSeqAlignment(new_alignments), output_file, "fasta")


# Remove divergent sequences
def remove_divergent(log, input_dir, output_dir, threshold, prefix, cores):
    temp_dir = tempfile.TemporaryDirectory()
    output_dir_div = make_dir(log, f"{output_dir}/{prefix}-03-divergent-removed")
    for file in glob.glob(os.path.join(input_dir, "*.fasta")):
        ucename = os.path.splitext(os.path.basename(file))[0]
        aln = AlignIO.read(file, "fasta")
        
        # Get the sequences as a 2D numpy array
        sequences = np.array([list(str(record.seq)) for record in aln])

        # Use multiprocessing to calculate similarities
        args_list = [(seq1, sequences) for seq1 in sequences]
        with Pool(cores) as p:
            similarities = p.map(calculate_similarity, args_list)

        # Filter sequences based on similarity threshold
        filtered_records = [record for record, similarity in zip(aln, similarities) if similarity >= threshold]

        # Write Report to log
        removed = [record for record, similarity in zip(aln, similarities) if similarity <= threshold]
        for record in removed:
            log.info(f"{record.id} from {ucename}")

        # Create a new MultipleSeqAlignment object with the filtered sequences
        filtered_alignment = MultipleSeqAlignment(filtered_records)

        # Write the new alignment to a file
        output_file = os.path.join(temp_dir.name, f"{ucename}.fasta")
        AlignIO.write(MultipleSeqAlignment(filtered_records), output_file, "fasta")

        
    # Use multiprocessing to calculate similarities
    args_list = [(output_dir_div, prefix, fasta_file) for fasta_file in glob.glob(os.path.join(temp_dir.name, "*.fasta"))]

    with Pool(cores) as p:
        p.map(remove_gap_columns, args_list)
    
    temp_dir.cleanup()

    return output_dir_div


# Filter Aligments #Taxa and bp Length
def filter_aligns(log, input_dir, output_dir, min_length, min_taxa, prefix):
    output_dir_filt = make_dir(log, f"{output_dir}/{prefix}-05-filtered")
    for file in glob.glob(os.path.join(input_dir, "*.fasta")):
        aln = AlignIO.read(file, "fasta")
        count = 0

        # Get min lenght of alignments
        if aln.get_alignment_length() >= min_length:
            length = True
        else:
            length = False

        # Get min number of taxa
        for taxon in aln:
            if (set(taxon.seq) == set("-")) or (set(taxon.seq) == set("?")):
                pass
            else:
                count += 1
        if count >= min_taxa:
            taxa = True
        else:
            taxa = False
    
        # Copy Files that met conditions
        if taxa and length:
            name = os.path.basename(file)
            shutil.copy(file, os.path.join(output_dir_filt, name))

    return output_dir_filt

# Get completeness matrices
def get_completeness(log, input_dir, output_dir, percent, taxa, prefix):
    output_dir_perc = make_dir(log, f"{output_dir}/{prefix}-{int(percent*100)}p")
    for file in glob.glob(os.path.join(input_dir, "*.fasta")):
        min_count = int(math.floor(percent * taxa))
        aln = AlignIO.read(file, "fasta")
        if len(aln) >= min_count:
            shutil.copyfile(file, os.path.join(output_dir_perc, os.path.basename(file)))
    return output_dir_perc

# Count Taxa in Alignments
def taxa_counter(input_dir, output_dir, output_format):
    output_filename = f"{output_dir}/summary-taxa.csv"
    count = Counter()
    for f in glob.glob(os.path.join(input_dir, f'*{output_format}')):
        aln = AlignIO.read(f, output_format)
        for seq in aln:
            if len(set(str(seq.seq))) > 1:
                count[seq.id] += 1
    with open(output_filename, "w") as outfile:
        outfile.write("taxon,count\n")
        for taxon, cnt in count.items():
            outfile.write("{},{}\n".format(taxon, cnt))

# Re-Align
def main():
    args = get_args()
    log = setup_logging(args)
    output_dir = make_dir(log, f"{args.prefix}-Alignments")
        
    # Align sequences
    align_command = (f"phyluce_align_seqcap_align "
        f"--input {args.monolithic_file} "
        f"--output {output_dir}/mafft-tmp "
        f"--taxa {args.taxa} --aligner mafft --cores {args.cores} "
        f"--output-format fasta "
        f"--incomplete-matrix "
        f"--ambiguous --no-trim"
        )
    text = f"Aligning sequences with MAFFT"
    log.info(text.center(65, "+"))
    subprocess.run(align_command, shell=True, check=True)

    tmp_alignments = f"{output_dir}/mafft-tmp"

    # Remove locus name
    log.info(f"Renaming alignments")
    alignments = batch_rename(log, tmp_alignments, output_dir, args.prefix, args.cores)
    shutil.rmtree(tmp_alignments)
    # Print a blank line for log formatting only
    print("")

    if args.trimAL is True:
        text = f"Trimming alignments with trimAL"
        log.info(text.center(55, "+"))
        alignments = trim_aligns(log, trimal_worker, alignments, output_dir, 
            args.t_method, args.b1, args.b2, args.b3, args.b4, 
            args.prefix, args.cores)

    if args.gblocks is True: # I Need to handle gblocks error here for not working on Mac for example
        text = f"Trimming alignments with gblocks"
        log.info(text.center(55, "+"))
        alignments = trim_aligns(log, gblocks_worker, alignments, output_dir, 
            args.t_method, args.b1, args.b2, args.b3, args.b4, 
            args.prefix, args.cores)

    if args.remove_div is True:
        text = f"Removing Divergent Sequences"
        log.info(text.center(55, "+"))
        alignments = remove_divergent(log, alignments, output_dir, args.divergent, args.prefix, args.cores)
        
    if args.remove_short is True:
        text = f"Removing Short Sequences"
        log.info(text.center(55, "+"))
        alignments = fasta_drop(log, alignments, output_dir, args.short_cutoff, args.prefix)

    if args.filter_alignments is True:
        text = f"Filtering Aligments"
        log.info(text.center(55, "+"))
        alignments = filter_aligns(log, alignments, output_dir, args.min_length, args.min_taxa, args.prefix)

    if args.get_completeness is True:
        text = f"Getting {int(args.percent*100)}% Completeness Alignments"
        log.info(text.center(55, "+"))
        alignments = get_completeness(log, alignments, output_dir, args.percent, args.taxa, args.prefix)

    if args.out_format is "fasta":
        text = f"fasta aligns written to {alignments}"
        log.info(text)
    else:
        alignments = convert_alignments(log, alignments, "fasta", output_dir, args.out_format, args.prefix)
        text = f"{args.out_format} aligns written to {alignments}"
        log.info(text)

    if args.taxa_count is True:
        text = f"Printing Taxa Count in Alignments to summary-taxa.csv"
        log.info(text)
        taxa_counter(alignments, output_dir, args.out_format)

    text = f" Completed FUSe "
    log.info(text.center(55, "="))

if __name__ == "__main__":
    main()