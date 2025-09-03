#!/usr/bin/env python3
"""
Script to harmonize sample IDs between VCF and phenotype files.
Removes '0_' prefix from sample IDs in VCF file to match phenotype format.
"""

import gzip
import sys
import os
import subprocess

def harmonize_vcf_sample_ids(input_vcf, output_vcf):
    """
    Harmonize VCF sample IDs by removing '0_' prefix to match phenotype format.
    
    Args:
        input_vcf (str): Path to input VCF file (.vcf.gz)
        output_vcf (str): Path to output harmonized VCF file (.vcf.gz)
    """
    
    print(f"Processing {input_vcf}...")
    
    # Use tabix to process the VCF file efficiently
    # First, get the header to see the sample IDs
    header_cmd = f"zcat {input_vcf} | grep '^#'"
    
    try:
        # Process header lines
        with gzip.open(output_vcf, 'wt') as outfile:
            with gzip.open(input_vcf, 'rt') as infile:
                for line in infile:
                    if line.startswith('#'):
                        if line.startswith('#CHROM'):
                            # This is the sample header line
                            parts = line.strip().split('\t')
                            if len(parts) > 9:  # Has sample columns
                                # Keep first 9 columns as-is (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
                                new_parts = parts[:9]
                                # Remove '0_' prefix from all sample IDs
                                for sample_id in parts[9:]:
                                    if sample_id.startswith('0_'):
                                        new_parts.append(sample_id[2:])  # Remove '0_' prefix
                                    else:
                                        new_parts.append(sample_id)
                                outfile.write('\t'.join(new_parts) + '\n')
                            else:
                                outfile.write(line)
                        else:
                            # Other header lines - write as-is
                            outfile.write(line)
                    else:
                        # Data lines - write as-is
                        outfile.write(line)
        
        print(f"Harmonized VCF file saved to {output_vcf}")
        print("Sample IDs now match phenotype file format (removed '0_' prefix)")
        
        # Create tabix index for the new file
        print("Creating tabix index...")
        subprocess.run(['tabix', '-p', 'vcf', output_vcf], check=True)
        print("Tabix index created successfully!")
        
    except Exception as e:
        print(f"Error processing file: {e}")
        # Clean up partial output file
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        raise

def main():
    if len(sys.argv) != 3:
        print("Usage: python harmonize_vcf_sample_ids.py <input_vcf_file> <output_harmonized_vcf_file>")
        print("Example: python harmonize_vcf_sample_ids.py GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01_harmonized.vcf.gz")
        sys.exit(1)
    
    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    
    if not os.path.exists(input_vcf):
        print(f"Error: Input VCF file {input_vcf} not found!")
        sys.exit(1)
    
    # Check if tabix is available
    try:
        subprocess.run(['tabix', '--help'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: tabix is not available. Please install tabix or ensure it's in your PATH.")
        sys.exit(1)
    
    try:
        harmonize_vcf_sample_ids(input_vcf, output_vcf)
        print("Successfully harmonized VCF file sample IDs!")
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 