import pysam


def get_sample_name_from_bam(bam: str) -> str:
    """
    Get the sample name available at the bam header
    """

    b = pysam.AlignmentFile(bam, "rb")
    header_dict = b.header.to_dict()

    sample_name = ""
    if "RG" in header_dict:
        if "SM" in header_dict["RG"][0]:
            sample_name = header_dict["RG"][0]["SM"]
    return sample_name

def reverse_complement(dna_sequence):
    # Define the complement mapping
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    # Generate the complement sequence
    complement_sequence = "".join(complement[base] for base in dna_sequence)

    # Reverse the sequence and return
    return complement_sequence[::-1]