import os
import sys
import pysam
import argparse
import pyximport
import re
import edlib
from collections import defaultdict
pyximport.install()
from src.assembler import DeBruijnAssembler
from src.bed import BedRecord
from src.vcf import VCFwriter, Record, VCFreader

def scan_complex_indels(active_regions: list, bam: str):
    '''
    '''

    indel_calls = list()
    samfile = pysam.AlignmentFile(bam, "rb")
    for region in active_regions:
        msg = f" INFO: scanning region {region}"
        print(msg)
        # spanning_r = get_spanning_coordinates_from_alns(bam, region.chr,
        #     region.start, region.end)
        list_variants = list()
        variants = dict()
        other_reads = list()
        for read in samfile.fetch(region.chr, region.start, region.end):
            if read.mapping_quality < 10:
                continue
            cx_indel = refine_cigar_operators(read)
            if cx_indel:
                if cx_indel['VAR_ID'] not in variants:
                    variants[cx_indel['VAR_ID']] = cx_indel
                    variants[cx_indel['VAR_ID']]['SEQ'].append(read)
                else:
                    variants[cx_indel['VAR_ID']]['SEQ'].append(read)
            if re.search(r"[0-9]+[S|M][0-9]+[S|M]",read.cigarstring):
                other_reads.append(read)
        for var in variants:
            read_list = variants[var]['SEQ'] + other_reads
            seq_list = list()
            read_list = variants[var]['SEQ'] + other_reads
            for read in read_list:
                seq_list.append(read.query_sequence)
            dbg = DeBruijnAssembler(seq_list, 21)
            contig_list = dbg.eulerian_walk()
            # print(contig_list)
            # if variants[var]['CHROM'] == "chr7":
            #     print(variants[var]['CHROM'] + " " + contig)
            if not contig_list:
                continue
            aln_stats = align_reads_to_contig(read_list, contig_list)
            depth_info = samfile.count_coverage(variants[var]['CHROM'],
                int(variants[var]['POS'])-1,
                int(variants[var]['POS']),
                quality_threshold = 0)
            depth = 0
            for base in depth_info:
                depth += base[0]
            if depth == 0:
                AF = 0
            else:
                AF = round(aln_stats['supporting_reads']/depth, 4)
            AC = aln_stats['supporting_reads']
            refc = depth-AC
            call_dict = {
                "CHROM": variants[var]['CHROM'],
                "POS": variants[var]['POS'],
                "REF": variants[var]['REF'],
                "ID": ".",
                "ALT": variants[var]['ALT'],
                "QUAL": ".",
                "FILTER": ".",
                "INFO": {
                    "SOURCE": "RivIndel",
                    "AC": AC,
                    "DP": depth,
                    "AF": AF
                },
                "FORMAT": "GT:AD:AF:DP:F1R2:F2R1:SB",
                "SAMPLE": f"0/1:{AC},{refc}:{AF}:{depth}:{aln_stats['fwd_reads']}:{aln_stats['rev_reads']}:."
            }
            indel_calls.append(call_dict)
    return indel_calls

def align_reads_to_contig(reads, contigs) -> dict:
    '''
    '''
    aln_dict = {
        'supporting_reads' : 0,
        'fwd_reads': 0,
        'rev_reads': 0
    }
    for contig in contigs:
        for read in reads:
            aln_stats = edlib.align(read.query_sequence, contig, mode = "HW", task = "path")
            if int(aln_stats['editDistance']) < 4:
                if read.is_forward:
                    aln_dict['fwd_reads'] += 1
                else:
                    aln_dict['rev_reads'] += 1
                aln_dict['supporting_reads'] += 1
    return aln_dict

def refine_cigar_operators(read: str) -> dict:
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
    var_dict = dict()
    if not read.cigartuples:
        return var_dict
    if len(read.cigartuples) > 1:
        # print(read)
        read_seq = read.query_sequence
        ref_seq = read.get_reference_sequence()
        op_str = ""
        for op in read.cigartuples:
            op_type = operation_dict[op[0]]
            op_span = op[1]
            if op_type == "H":
                continue
            for i in range(0, op_span):
                op_str+=op_type

        valid_ops = ["D", "I", "M"]
        num = 0
        for op in valid_ops:
            if not op in op_str:
                pass
            else:
                num+=1
        if num < 2:
            return var_dict

        idx = 0
        jdx = 0
        read_cigar_list = list()
        ref_cigar_list = list()
        for op in op_str:
            if op == "H":
                continue
            if op == "M":
                read_cigar_list.append(read_seq[idx].upper())
                ref_cigar_list.append(ref_seq[jdx].upper())
                idx+=1
                jdx+=1
            if op == "D":
                read_cigar_list.append("*")
                ref_cigar_list.append(ref_seq[jdx].upper())
                jdx+=1
            if op == "I" or op == "S":
                read_cigar_list.append(read_seq[idx].upper())
                ref_cigar_list.append("*")
                idx+=1

        flag = False
        start = -1
        max_distance = 10
        end = -5
        for i in range(0, len(op_str)):

            ref_ntd = ref_cigar_list[i].upper()
            read_ntd = read_cigar_list[i].upper()
            if flag is False:
                if (ref_ntd != read_ntd and op_str[i] == "M") or (op_str[i] == "D" or op_str[i] == "I"):
                    flag = True
                    start = i
                    end = i
            if flag is True:
                if (ref_ntd != read_ntd and op_str[i] == "M") or (op_str[i] == "D" or op_str[i] == "I"):
                    if i-end <= max_distance:
                        if i > end:
                            end = i
                    else:
                        start = i
                        end = i
        if start > 0:
            if (end - start) >= 2:
                ref = ''.join(filter(lambda char: char != '*', ref_cigar_list[start:end+1]))
                alt = ''.join(filter(lambda char: char != '*', read_cigar_list[start:end+1]))
                var_id = f"{read.pos+start+1}-{ref}-{alt}"
                var_dict = {
                    "VAR_ID": var_id,
                    "CHROM": f"chr{str(read.reference_id)}",
                    "POS": read.pos+start+1,
                    "REF": ref,
                    "ALT": alt,
                    "SEQ": list()
                }
    return var_dict

def get_spanning_coordinates_from_alns(bam: str, chr: str, start: int,
    end: int)-> BedRecord:
    '''
        Given a BED record, loop throught overlapping reads and get the minimum
        and maximum genomic coordinates.
    '''
    max_start = 10e10
    max_end = 0
    samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile.fetch(chr, start, end):
        if read.pos < max_start:
            max_start = read.pos
        if not read.is_mapped:
            continue
        if read.mapping_quality < 10:
            continue
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
    vcf_out = args.vcf_out
    min_reads = args.min_reads

    active_regions = parse_active_regions(bed)

    vcf = VCFwriter(vcf_out, sample_name="TEST")
    header = vcf.header
    vcf.write(header)
    format = "GT:AD:AF:DP:F1R2:F2R1:SB"
    indel_calls = scan_complex_indels(active_regions, bam)
    for indel in indel_calls:

        if indel['INFO']['AC'] < min_reads:
            continue

        info_list = list()
        for field in indel['INFO']:
            info_list.append(f"{field}={str(indel['INFO'][field])}")
        info_str = ';'.join(info_list)

        out_list = [indel['CHROM'], str(indel['POS']),
            indel['ID'], indel['REF'],indel['ALT'],
            indel['QUAL'], indel['FILTER'], info_str, format, indel['SAMPLE']]
        out_str = '\t'.join(out_list)

        vcf.write(out_str)
    vcf.close()


def parse_arguments():
    '''
    '''
    parser = argparse.ArgumentParser(description='Description: Detect complex Indels given a list of regions')
    parser.add_argument( '--bam', dest='bam_file', type = str,
        help="Input BAM file", required=True)
    parser.add_argument( '--bed', dest='bed_file', type = str,
        help="BED regions to search", required=True)
    parser.add_argument( '--vcf', dest='vcf_out', type = str,
        help="Output vcf", required=True)
    parser.add_argument( '--min_reads', dest='min_reads', type = int,
        help="Minimum read support", default=2)


    args = parser.parse_args()
    return args

if __name__ == '__main__':

    main()
