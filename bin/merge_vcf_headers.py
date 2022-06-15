#!/usr/bin/env python

import argparse
import re
import logging
import sys

logger = logging.getLogger()

if __name__ == "__main__":
    # Define and parse the arguments
    parser = argparse.ArgumentParser(description="A script to add the pedigree and sample metadata of an empty VCF to another VCF")
    parser.add_argument('vcf1', metavar='VCF1', type=str, help="The first VCF file")
    parser.add_argument('vcf2', metavar='VCF2', type=str, help="The second VCF file")
    parser.add_argument('output', metavar='OUTPUT', type=str, help="The file to output the merged VCF to")

    args = parser.parse_args()

    vcf1 = args.vcf1
    vcf2 = args.vcf2
    output = args.output

    # Open and read the input files, open the output file
    file_vcf1 = open(vcf1, "r")
    file_vcf2 = open(vcf2, "r")
    output_file = open(output, "w")

    read_vcf1 = file_vcf1.read()
    read_vcf2 = file_vcf2.read()

    # Some quick checks to see if the files are compatible
    pedigree_pattern = '##PEDIGREE=.*\s'
    sample_pattern = '##SAMPLE=.*\s'

    if re.search(pedigree_pattern, read_vcf1):
        pedigree_file = read_vcf1
        vcf = read_vcf2
    elif re.search(pedigree_pattern, read_vcf2):
        pedigree_file = read_vcf2
        vcf = read_vcf1
    else:
        logger.critical("No '##PEDIGREE' header found inside both files")
        sys.exit(1)

    header_pattern = '#CHROM.*\s'
    if re.findall(header_pattern, pedigree_file) != re.findall(header_pattern, vcf):
        logger.critical("The header line is not the same in both files")
        sys.exit(1)

    # Find the pedigree and sample lines
    pedigree = re.findall(pedigree_pattern, pedigree_file)
    sample = re.findall(sample_pattern, pedigree_file)

    # Write the new VCF file
    status = "before_info"
    info_pattern = '^##INFO.*$'

    for line in vcf.split("\n"):
        if status == "before_info" and re.search(info_pattern, line):
            status = "during_info"
        elif status == "during_info" and not re.search(info_pattern, line):
            output_file.writelines(sample)
            output_file.writelines(pedigree)
            status = "after_info"
        output_file.write(f'{line}\n')

    # Close all files
    file_vcf1.close()
    file_vcf2.close()
    output_file.close()
