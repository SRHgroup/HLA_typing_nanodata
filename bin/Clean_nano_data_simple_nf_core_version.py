#!/usr/bin/env python

import argparse
import csv
import regex


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict[base] for base in reverse_sequence)
    return reverse_complement_sequence

def find_primer_index_fast(sequence, primer, max_mismatches=2):
    primer_len = len(primer)

    for i in range(len(sequence) - primer_len + 1):
        mismatch_count = 0
        for j in range(primer_len):
            if sequence[i + j] != primer[j]:
                mismatch_count += 1
                if mismatch_count > max_mismatches:
                    break  # No need to check further if max mismatches exceeded

        if mismatch_count <= max_mismatches:
            return i

    return -1  # Primer not found


def fuzzyMapFind(subSeq, sequence, subs=2):
    pattern = "(?e)({}){{s<={}}}".format(subSeq, subs)
    match = regex.search(pattern, sequence)
    if match is None:
        return -1
    else:
        return match.start()

def findPrimerMatchesNested(forwardPrimers,sequence,quality,delta=500):
    reversePrimers  = [reverse_complement(fp) for fp in forwardPrimers]
    for fp in forwardPrimers:
        #fp_idx = find_primer_index_fast(sequence,fp)
        fp_idx = fuzzyMapFind(fp,sequence)
        for rp in reversePrimers:
            #rp_idx = find_primer_index_fast(sequence,rp)
            rp_idx = fuzzyMapFind(rp,sequence)
            if fp_idx >= 0 and rp_idx >= 0 and rp_idx-fp_idx>delta:
                trimmed_sequence = sequence[fp_idx:rp_idx+len(fp)]
                trimmed_quality = quality[fp_idx:rp_idx+len(fp)]
                return (trimmed_sequence,trimmed_quality)
                # Write the modified sequence entry to the output file
                return (trimmed_sequence,trimmed_quality)     
            elif fp_idx >= 0 and rp_idx < 0:
                trimmed_sequence = sequence[fp_idx:len(sequence)]
                trimmed_quality = quality[fp_idx:len(sequence)]
                # Write the modified sequence entry to the output file
                return (trimmed_sequence,trimmed_quality)                     
    return (None,None)



def trim_fastq(input_file, output_file, junkfile,countfile):
    locus_counts = {locus: {'matches': 0, 'no_matches': 0} for locus in primer_dict}
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile, open(junkfile, 'w') as junk:#, open(countfile, 'w') as count
        while True:
            # Read four lines for each sequence entry in the FASTQ file
            header = infile.readline().strip()
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()
            quality = infile.readline().strip()

            # Check if we have reached the end of the file
            if not header:
                break
         #   outputGather = []
            
            for loci,forwardPrimers in primer_dict.items():
             #   print(loci)
             #   print(forwardPrimers)
                trimmedOutput = findPrimerMatchesNested(forwardPrimers,sequence,quality,delta=300)
                if not trimmedOutput[0] is None:
                    locus_counts[loci]['matches'] += 1
                    outfile.write(f"{header}\n{trimmedOutput[0]}\n{plus_line}\n{trimmedOutput[1]}\n")
                    #outputGather.append(trimmedOutput)
                else:
                    locus_counts[loci]['no_matches'] += 1
                    junk.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")
                   # outfile.write(f"{header}\n{sequence}\n{plus_line}\n{quality}\n")

    # Print locus counts

def read_primer_csv(file_path):
    primer_dict = {}
    with open(file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter=';')
        for row in csv_reader:
            primer_name = row['PrimerName']
            primer_prefix = primer_name.split('-')[0] 
            primer_sequence = row['PrimerSequence'].upper()
            if primer_prefix in primer_dict:
                primer_dict[primer_prefix].append(primer_sequence)
            else:
                primer_dict[primer_prefix] = [primer_sequence]
    return primer_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean nonopore files based on primers")
    parser.add_argument("--file", "-i", type=str, required=True, help="Input file path")
    parser.add_argument("--primer","-p", required=True,help="CSV file containing primer information")
    parser.add_argument("--outfile", "-o", type=str, required=True, help="Output file name")
    parser.add_argument("--junkfile", "-j", type=str, required=True, help="junk file name")
    parser.add_argument("--countfile", "-c", type=str, required=True, help="lucus count file name")
    args = parser.parse_args()



# Warnings for input files
if not args.file:
    warnings.warn("Input file not provided. Please specify the --file or -i argument.")
if not args.primer:
    warnings.warn("Input primer not provided. Please specify the --primer or -p argument.")
if not args.outfile:
    warnings.warn("Output file not provided. Please specify the --outfile or -o argument.")
if not args.junkfile:
    warnings.warn("Junk file not provided. Please specify the --junkfile or -j argument.")
if not args.countfile:
    warnings.warn("count file not provided. Please specify the --countfile or -c argument.")

#aign parsed values 
input_fastq = args.file
primer_csv = args.primer
output_fastq = args.outfile
junk = args.junkfile
count = args.countfile

primer_dict = read_primer_csv(primer_csv)
trim_fastq(input_fastq, output_fastq, junk,count)





# how to run 
#python3 bin/Clean_Nano_Data_all_primers.py -i barcode12/PAS44137_pass_barcode12_bdac54ca_131b8794_0.fastq -o barcode12/out_test_barcode12.fastq -j barcode12/junk_barcode12.fastq
