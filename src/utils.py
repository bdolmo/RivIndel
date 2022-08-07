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
