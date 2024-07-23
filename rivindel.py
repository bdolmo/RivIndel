import argparse
import subprocess
import os

script_dir = os.path.dirname(os.path.realpath(__file__))


def run_rivindel(args):
    # Determine the directory where the script is located
    
    # Construct the path to the binary relative to the script's directory
    binary_path = os.path.join(script_dir, 'src', 'rivindel')
    
    # Construct the command
    command = [
        binary_path,
        '--bam', args.bam,
        '--threads', str(args.threads),
        '--ref', args.ref,
        '--vcf', args.vcf,
        '--bed', args.bed
    ]
    
    # Execute the command
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Command output:", result.stdout.decode())
        print("Command error:", result.stderr.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        print("Command output:", e.stdout.decode())
        print("Command error:", e.stderr.decode())

def normalize_vcf(input_vcf, output_vcf, reference_fasta):
    """ """
    bcftools_path = os.path.join(script_dir, "bcftools-1.20", "bcftools")

    if not os.path.isfile(bcftools_path):
        raise FileNotFoundError(f"bcftools executable not found at {bcftools_path}")
    if not os.access(bcftools_path, os.X_OK):
        raise PermissionError(f"bcftools executable at {bcftools_path} is not executable")

    command = [
        bcftools_path, 'norm',
        '-f', reference_fasta,
        '-o', output_vcf,
        '-O', 'v',
        input_vcf
    ]
    
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("Normalization output:", result.stdout.decode())
        print("Normalization error:", result.stderr.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error normalizing VCF: {e}")
        print("Normalization output:", e.stdout.decode())
        print("Normalization error:", e.stderr.decode())

def is_file_empty(file_path):
    return os.path.exists(file_path) and os.path.getsize(file_path) == 0

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run RivIndel binary with specified arguments.")
    
    # Define arguments
    parser.add_argument('--bam', required=True, help='Path to the BAM file')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--ref', required=True, help='Path to the reference fasta file')
    parser.add_argument('--vcf', required=True, help='Path to the output VCF file')
    parser.add_argument('--bed', help='Path to the BED file')
    parser.add_argument('--no-norm-vcf', action='store_false', dest='norm_vcf', help='Skip normalization of the VCF file')
    parser.add_argument('--force', action='store_true', help='Force execution of the RivIndel binary')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if the output VCF exists and is not empty, or if force is specified
    if args.force or not os.path.exists(args.vcf) or is_file_empty(args.vcf):
        pass
        # run_rivindel(args)
    else:
        print(f"Skipping RivIndel execution as the output VCF file {args.vcf} aleady exists")
    
    # Normalize the VCF if --no-norm-vcf is not provided
    if args.norm_vcf:
        normalized_vcf_path = args.vcf.replace('.vcf', '.normalized.vcf')
        normalize_vcf(args.vcf, normalized_vcf_path, args.ref)

if __name__ == "__main__":
    main()
