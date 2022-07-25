import os
import sys
import pysam
import argparse


class BedRecord():
    '''
        Class to store BED records.
    '''
    def __init__(self, chr: str, start: int, end: int):
        self.chr = chr
        self.start = start
        self.end = end
        if start > end:
            msg = f" ERROR: start {start} cannot be greater than end {end}"
            raise ValueError(msg)
    def __repr__(self):
        return f"{self.chr}\t{self.start}\t{self.end}"


def scan_indels(active_regions: list, bam: str, genome_fasta: str):
    '''
    '''

    samfile = pysam.AlignmentFile(bam, "rb")
    for region in active_regions:

        spanning_r = get_spanning_coordinates_from_alns(bam, region.chr,
            region.start, region.end)

        for read in samfile.fetch(region.chr, region.start, region.end):
            # ref_sequence = read.get_reference_sequence()
            refine_cigar_operators(read)

def refine_cigar_operators(read):
    '''
    '''
    operation_dict = {
        0: "M",
        1: "I",
        2: "D",
        3: "N",
        4: "S",
        5: "H",
        6: "P",
        7: "=",
        8: "X",
        9: "B"
    }

    if len(read.cigartuples) > 1:
        read_seq = read.query_sequence
        ref_seq = read.get_reference_sequence()
        # print(ref_seq)
        # print(read_seq)
        op_str = ""
        for op in read.cigartuples:
            op_type = operation_dict[op[0]]
            op_span = op[1]
            for i in range(0, op_span):
                op_str+=op_type
        idx = 0
        jdx = 0
        read_str = ""
        ref_str = ""
        for op in op_str:
            if op == "M":
                read_str+=read_seq[idx]
                ref_str+=ref_seq[jdx]
                idx+=1
                jdx+=1
            if op == "D":
                read_str+= "-"
                ref_str+=ref_seq[jdx]
                jdx+=1

            if op == "I":
                ref_str+= "-"
                read_str+=read_seq[idx]
                idx+=1
        print(ref_str)
        print(op_str)
        print(read_str)
        print("\n")



        # sys.exit()
        # print(op_str)
        # print(read_str)
        # print("\n")


def get_spanning_coordinates_from_alns(bam: str, chr: str, start: int,
    end: int)-> BedRecord:
    '''
        Give a BED record, traverse all overlapping reads and retrieve the
        entire overlapping region
    '''
    max_start = 10e10
    max_end = 0
    samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile.fetch(chr, start, end):
        if read.pos < max_start:
            max_start = read.pos
        read_end = read.pos + read.reference_length
        if read_end > max_end:
            max_end = read_end
    bed_r = BedRecord(chr, max_start, max_end)
    return bed_r

def parse_active_regions(bed: str) -> list():
    '''
        Read a givem BED and return a list of regions to operate.
        An evolution of this function could be to parse the entire BAM file
        to locate candidate complex indels.
    '''

    active_regions = []
    with open(bed) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            tmp = line.split("\t")
            if len(tmp) < 3:
                msg = f" ERROR: line {line} cannot be parsed into at least three items"
                raise ValueError(msg)
            chr = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])
            bed_r = BedRecord(chr, start, end)
            active_regions.append(bed_r)
    f.close()
    return active_regions

def main():
    '''
    '''
    args = parse_arguments()
    bam = args.bam_file
    bed = args.bed_file
    genome_fasta = args.genome_fasta

    active_regions = parse_active_regions(bed)
    scan_indels(active_regions, bam, genome_fasta)

def parse_arguments():
    '''
    '''
    parser = argparse.ArgumentParser(description='Description: Detect complex Indels given a list of regions')
    parser.add_argument( '--bam', dest='bam_file', type = str,
        help="Input BAM file", required=True)
    parser.add_argument( '--bed', dest='bed_file', type = str,
        help="Regions to operate", required=True)
    parser.add_argument( '--fasta', dest='genome_fasta', type = str,
        help="Genome file in FASTA format", required=True)

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
