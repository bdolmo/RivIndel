# import pytest
import pyximport

pyximport.install()
from src.rivindel.scan import parse_active_regions, scan_complex_indels
from src.rivindel.report import write_variants


def test_positive_case():
    """
    Test whether a known variant can be detected from a test bam file
    """

    variant = "chr7\t55242467\t.\tAATTAAGAGAAG\tAGC"

    bam = "test/test.bam"
    bed = "test/test.bed"
    vcf_out = "test/test.vcf"
    min_reads = 4

    active_regions = parse_active_regions(bed)
    indel_calls = scan_complex_indels(active_regions, bam)
    write_variants(indel_calls, vcf_out, bam, min_reads)

    with open(vcf_out) as fh:
        results = fh.read()

    assert variant in results
