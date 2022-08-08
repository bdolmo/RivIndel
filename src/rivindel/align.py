import edlib


def align_reads_to_contigs(reads, contigs) -> dict:
    """ """
    aln_dict = {"supporting_reads": 0, "fwd_reads": 0, "rev_reads": 0}
    for contig in contigs:
        for read in reads:
            aln_stats = edlib.align(read.query_sequence, contig, mode="HW", task="path")
            if int(aln_stats["editDistance"]) < 4:
                if read.is_forward:
                    aln_dict["fwd_reads"] += 1
                else:
                    aln_dict["rev_reads"] += 1
                aln_dict["supporting_reads"] += 1
    return aln_dict
