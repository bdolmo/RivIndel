import argparse
import pyximport
pyximport.install()

from src.rivindel.scan import parse_active_regions, scan_complex_indels
from src.rivindel.report import write_variants


def main():
    """ """
    args = parse_arguments()
    bam = args.bam_file
    bed = args.bed_file
    vcf_out = args.vcf_out
    min_reads = args.min_reads

    active_regions = parse_active_regions(bed)
    indel_calls = scan_complex_indels(active_regions, bam)
    write_variants(indel_calls, vcf_out, bam, min_reads)


def parse_arguments():
    """ """
    parser = argparse.ArgumentParser(
        description="Description: Detect complex Indels given a list of regions"
    )
    parser.add_argument(
        "--bam", dest="bam_file", type=str, help="Input BAM file", required=True
    )
    parser.add_argument(
        "--bed", dest="bed_file", type=str, help="BED regions to search", required=True
    )
    parser.add_argument(
        "--vcf", dest="vcf_out", type=str, help="Output vcf", required=True
    )
    parser.add_argument(
        "--min_reads",
        dest="min_reads",
        type=int,
        help="Minimum read support",
        default=2,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    main()
