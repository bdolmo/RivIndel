from .utils import get_sample_name_from_bam
from .vcf import VCFwriter


def write_variants(indel_calls: dict, vcf_out: str, bam: str, min_reads: int) -> None:
    """ """
    sample_name = get_sample_name_from_bam(bam)

    vcf = VCFwriter(vcf_out, sample_name=sample_name)
    header = vcf.header
    vcf.write(header)
    fmt = "GT:AD:AF:DP:F1R2:F2R1:SB"

    for indel in indel_calls:
        # skip low allele count variants
        if indel["INFO"]["AC"] < min_reads:
            continue
        info_list = list()
        for field in indel["INFO"]:
            info_list.append(f"{field}={str(indel['INFO'][field])}")
        info_str = ";".join(info_list)

        out_list = [
            indel["CHROM"],
            str(indel["POS"]),
            indel["ID"],
            indel["REF"],
            indel["ALT"],
            indel["QUAL"],
            indel["FILTER"],
            info_str,
            fmt,
            str(indel["SAMPLE"]),
        ]
        out_str = "\t".join(out_list)

        vcf.write(out_str)
    vcf.close()
