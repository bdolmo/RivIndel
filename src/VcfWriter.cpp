#include "VcfWriter.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>


VcfWriter::VcfWriter(const std::string& vcfFilename, const std::string& bamFilename)
    : samFilePtr(nullptr), bamHeader(nullptr) {
    vcfFile.open(vcfFilename);
    if (!vcfFile.is_open()) {
        throw std::runtime_error("Failed to open VCF file " + vcfFilename);
    }

    samFilePtr = sam_open(bamFilename.c_str(), "r");
    if (!samFilePtr) {
        throw std::runtime_error("Failed to open BAM file " + bamFilename);
    }

    bamHeader = sam_hdr_read(samFilePtr);
    if (!bamHeader) {
        sam_close(samFilePtr);
        throw std::runtime_error("Failed to read header from BAM file " + bamFilename);
    }

    // Default name 
    sampleName = "SAMPLE";

    std::string headerAsString = std::string(bamHeader->text, bamHeader->l_text);
    std::istringstream headerStream(headerAsString);
    std::string line;
    while (std::getline(headerStream, line)) {
        if (line.substr(0, 3) == "@RG") {
            size_t smPos = line.find("SM:");
            if (smPos != std::string::npos) {
                smPos += 3;
                size_t smEnd = line.find('\t', smPos);
                if (smEnd == std::string::npos) {
                    sampleName = line.substr(smPos);
                } 
                else {
                    sampleName = line.substr(smPos, smEnd - smPos);
                }
            }
        }
    }
}


VcfWriter::~VcfWriter() {
    if (vcfFile.is_open()) {
        vcfFile.close();
    }
    if (bamHeader) {
        bam_hdr_destroy(bamHeader);
    }
    if (samFilePtr) {
        sam_close(samFilePtr);
    }
}
void VcfWriter::writeHeader() {

    vcfFile << "##fileformat=VCFv4.3\n";
    vcfFile << "##source=RivIndel\n";

    if (bamHeader) {
        for (int i = 0; i < bamHeader->n_targets; ++i) {
            vcfFile << "##contig=<ID=" << bamHeader->target_name[i]
                    << ",length=" << bamHeader->target_len[i] << ">\n";
        }
    }
    vcfFile << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n";
    vcfFile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n";
    vcfFile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n";
    vcfFile << "##INFO=<ID=FWD,Number=1,Type=Integer,Description=\"Forward strand supporting reads\">"<< "\n";
    vcfFile << "##INFO=<ID=REV,Number=1,Type=Integer,Description=\"Reverse strand supporting reads\">"<< "\n";
    vcfFile << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias (0 to 1)\">"<< "\n";
    vcfFile << "##INFO=<ID=KDIV_CONTIG,Number=1,Type=Float,Description=\"k-mer diversity of the contig\">"<< "\n";
    vcfFile << "##INFO=<ID=KDIV_5Flank,Number=1,Type=Float,Description=\"k-mer diversity of the 5' Upstream sequence\">"<< "\n";   
    vcfFile << "##INFO=<ID=KDIV_3Flank,Number=1,Type=Float,Description=\"k-mer diversity of the 3' Downstream sequence\">"<< "\n";
    vcfFile << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">"<< "\n";   
    vcfFile << "##INFO=<ID=MBQ,Number=1,Type=Float,Description=\"Mean Base Quality\">"<< "\n";
    vcfFile << "##INFO=<ID=GC,Number=1,Type=Float,Description=\"GC-content of the contig\">"<< "\n";
    vcfFile << "##INFO=<ID=HOML,Number=1,Type=Integer,Description=\"Longest homopolymer run\">"<< "\n";   
    vcfFile << "##INFO=<ID=HOMN,Number=1,Type=Integer,Description=\"Number of homopolymers\">"<< "\n";   
    vcfFile << "##INFO=<ID=DNTDL,Number=1,Type=Integer,Description=\"Longest di-nucleotide run\">"<< "\n";   
    vcfFile << "##INFO=<ID=DNTDN,Number=1,Type=Integer,Description=\"Number of di-nucleotide runs\">"<< "\n";  
    vcfFile << "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Caller name\">"<< "\n";

    vcfFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n";
    vcfFile << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"<< "\n";
    vcfFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth of Coverage of the variant\">"<< "\n";
    vcfFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
}

void VcfWriter::writeVariant(const variant_t& var) {

    vcfFile << var.chr << "\t"
            << var.pos<< "\t"
            << ".\t"  // ID 
            << var.ref << "\t"
            << var.alt << "\t"
            << ".\t"  // QUAL
            << ".\t"
            << "AC=" << var.readSupport 
            << ";AF=" << var.alleleFrequency 
            << ";DP=" << var.totalDepth
            << ";FWD=" << var.plusStrand
            << ";REV=" << var.minusStrand
            << ";SB=" << var.strandBias
            << ";KDIV_CONTIG=" << var.kmerDiversity
            << ";KDIV_5Flank=" << var.flankingKmerDiversityUpstream
            << ";KDIV_3Flank=" << var.flankingKmerDiversityDownstream
            << ";MQ=" << var.mapQual
            << ";MBQ=" << var.meanBaseQuality
            << ";GC=" << var.gcContent
            << ";HOML=" << var.longestHomopolymerRun
            << ";HOMN=" << var.numberOfHomopolymerRuns
            << ";DNTDL=" << var.longestDintdRun
            << ";DNTDN=" << var.numberOfDintd
            << ";SOURCE=RivIndel"<< "\t"
            << "GT:AD:DP" << "\t"
            << "./." << ":" << var.alleleDepth << ":" << var.totalDepth
            << "\n"; 
}
