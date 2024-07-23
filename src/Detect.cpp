#include <iostream>
#include <string>
#include <htslib/sam.h>
#include "BamRecord.h"
#include "BamReader.h"
#include <vector>
#include <optional>
#include "RefFasta.h"
#include "sw.cpp"
#include "FermiAssembler.cpp"
#include "GreedyAssembler.cpp"
#include "VcfWriter.h"
#include <utility>
#include "variant_t.h"
#include <cctype>
#include <numeric>
#include <algorithm>
#include "../edlib/include/edlib.h"
#include <unordered_set>
#include <utility>
#include <omp.h>
#include <cmath> 
#include <mutex>
#include "ThreadSafeQueue.h"
// #include "olc.cpp"
#include "Assembler.h"
// #include "vcake.cpp"
#include <unordered_set>


double ComputeStrandBias(int forwardStrandCount, int reverseStrandCount) {
    if (forwardStrandCount + reverseStrandCount == 0) {
        return 0.0;
    }
    
    return static_cast<double>(std::fabs(forwardStrandCount - reverseStrandCount)) / (forwardStrandCount + reverseStrandCount);
}

struct MutationSegment {
    int start;
    int end;
};


std::pair<int, int> calculateDinucleotideMetrics(const std::string& sequence, int minLength = 4) {
    int longestRun = 0;
    int currentRunLength = 1;
    int numberOfRuns = 0;

    if (sequence.length() < 2) {
        return {0, 0}; // Return if the sequence is too short for dinucleotide runs
    }

    // Helper function to check if a dinucleotide is not a homopolymer
    auto isDifferentDinucleotide = [](const std::string& dinuc) {
        return dinuc.length() == 2 && dinuc[0] != dinuc[1];
    };

    for (size_t i = 2; i < sequence.length(); i += 2) {
        std::string currentDinuc = sequence.substr(i, 2);
        std::string previousDinuc = sequence.substr(i - 2, 2);

        if (isDifferentDinucleotide(currentDinuc) && currentDinuc == previousDinuc) {
            currentRunLength++;
        } else {
            if (currentRunLength >= minLength) {
                numberOfRuns++;
                longestRun = std::max(longestRun, currentRunLength);
            }
            currentRunLength = isDifferentDinucleotide(currentDinuc) ? 1 : 0;
        }
    }

    // Check the last run
    if (currentRunLength >= minLength) {
        numberOfRuns++;
        longestRun = std::max(longestRun, currentRunLength);
    }

    return {longestRun, numberOfRuns};
}



// Function to calculate homopolymer run metrics with a minimum length of 8
std::pair<int, int> calculateHomopolymerMetrics(const std::string& sequence, int minLength = 8) {
    int longestRun = 0;
    int currentRunLength = 1;
    int numberOfRuns = 0;

    for (size_t i = 1; i < sequence.length(); ++i) {
        if (sequence[i] == sequence[i - 1]) {
            currentRunLength++;
        } else {
            if (currentRunLength >= minLength) {
                numberOfRuns++;
                longestRun = std::max(longestRun, currentRunLength);
            }
            currentRunLength = 1;
        }
    }

    // Check the last run
    if (currentRunLength >= minLength) {
        numberOfRuns++;
        longestRun = std::max(longestRun, currentRunLength);
    }

    return {longestRun, numberOfRuns};
}

float calculateGCContent(const std::string& sequence) {
    int gcCount = 0;
    for (char base : sequence) {
        if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
            gcCount++;
        }
    }
    float gcContent = static_cast<float>(gcCount) / sequence.length();
    return gcContent;
}


// Function to calculate k-mer diversity
float calculateKmerDiversity(const std::string& sequence, int k=3) {
    std::unordered_set<std::string> uniqueKmers;
    int n = sequence.size();

    // Extract all k-mers and add to the set
    for (int i = 0; i <= n - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        uniqueKmers.insert(kmer);
    }

    // Calculate diversity as the number of unique k-mers divided by the total possible k-mers
    int totalKmers = n - k + 1;
    if (totalKmers <= 0) return 0.0f; // Handle edge case where sequence length < k
    float diversity = static_cast<float>(uniqueKmers.size()) / totalKmers;

    return diversity;
}

// Function to extract flanking sequences and calculate k-mer diversity
std::pair<float, float> getFlankingKmerDiversity(const std::string& chr, int pos, const std::string& refFasta, int flankSize = 50) {
    RefFasta ref(refFasta);
    
    // Extract upstream and downstream sequences from the reference genome
    std::string upstream = ref.fetchSequence(chr, pos - flankSize, pos - 1);
    std::string downstream = ref.fetchSequence(chr, pos + 1, pos + flankSize);

    // Calculate k-mer diversity for the flanking sequences
    float upstreamDiversity = calculateKmerDiversity(upstream, 3);
    float downstreamDiversity = calculateKmerDiversity(downstream, 3);

    return {upstreamDiversity, downstreamDiversity};
}



std::pair<int, std::string> alignReadToContig(const std::string& read, const std::string& contig) {
    EdlibAlignResult result = edlibAlign(read.c_str(), read.size(), contig.c_str(), contig.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

    int editDistance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
    std::string cigar;

    if (result.status == EDLIB_STATUS_OK && result.alignment) {
        char* cigarC = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        cigar = std::string(cigarC);
        free(cigarC);
    }

    edlibFreeAlignResult(result);
    return {editDistance, cigar};
}


// Function to find the complement of a nucleotide
char complement(char nucleotide) {
    switch (nucleotide) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default: return nucleotide; // For non-ATCG characters, return as is
    }
}

std::string reverseComplement(const std::string& dna) {
    std::string revComp(dna.rbegin(), dna.rend()); // Reverse the string
    std::transform(revComp.begin(), revComp.end(), revComp.begin(), complement);
    return revComp;
}

struct contig_t{
    std::string seq;
    std::string chr;
    int start;
    int end;
    int readSupport;
    int plusStrand;
    int minusStrand;
};

float calculateMean(const std::vector<int>& nums) {
    if (nums.empty()) {
        return 0;
    }
    int sum = std::accumulate(nums.begin(), nums.end(), 0);
    float mean = static_cast<float>(sum) / nums.size();
    return mean;
}

int ExtractChromosomeNumber(const std::string& chr) {
    size_t i = 0;
    while (i < chr.length() && !std::isdigit(chr[i])) ++i;
    size_t start = i;
    while (i < chr.length() && std::isdigit(chr[i])) ++i;
    if (start < i) {
        return std::stoi(chr.substr(start, i - start));
    }
    return 0;
}

bool CompareVariants(const variant_t& a, const variant_t& b) {
    int numA = ExtractChromosomeNumber(a.chr);
    int numB = ExtractChromosomeNumber(b.chr);

    if (numA != 0 && numB != 0) {
        if (numA != numB) return numA < numB;
    } 
    else {
        if (a.chr != b.chr) return a.chr < b.chr;
    }

    return a.pos < b.pos;
}

void PopulateVCF(const std::vector<variant_t>& variants, const std::string& bamFile, const std::string& vcfFilename) {
    std::vector<variant_t> allVariants = variants;

    // Sort all variants by chromosome and position
    std::sort(allVariants.begin(), allVariants.end(), CompareVariants);

    // Create VcfWriter to write the VCF file
    VcfWriter writer(vcfFilename, bamFile);
    writer.writeHeader();
    writer.writeVariants(allVariants);
}


// void PopulateVCF(const std::vector<variant_t>& variants, const std::string& bamFile, const std::string& vcfFilename) {
//     std::vector<variant_t> allVariants = variants;

//     std::sort(allVariants.begin(), allVariants.end(), CompareVariants);

//     VcfWriter writer(vcfFilename, bamFile);
//     writer.writeHeader();
//     for (const auto& var : allVariants) {
//         writer.writeVariant(var);
//     }
// }

std::vector<BamRecord> GetSoftClippedReads(BamReader& reader) {
    std::vector<BamRecord> softClippedReads;
    BamRecord record;

    while (reader.GetNextRecord(record)) {
        std::string cigar = record.cigarString();
        size_t cigarLength = cigar.length();
        if (cigarLength > 1) { 
            int leadingSoftClips = 0;
            int trailingSoftClips = 0;
            bool leading = true;

            for (size_t i = 0; i < cigarLength; ++i) {
                char c = cigar[i];
                if (isdigit(c)) { 
                    (leading ? leadingSoftClips : trailingSoftClips) = (leading ? leadingSoftClips : trailingSoftClips) * 10 + (c - '0');
                } 
                else if (c == 'S') {
                    if (i > 0 && !leading) trailingSoftClips = 0;
                    leading = false;
                } 
                else if (c == 'M') {
                    if (!leading) break;
                    leading = false;
                } 
                else {
                    leadingSoftClips = 0;
                    trailingSoftClips = 0;
                    break;
                }
            }

            if (leadingSoftClips >= 10 || trailingSoftClips >= 10) {
                softClippedReads.push_back(record);
            }
        }
    }
    return softClippedReads;
}


int CalculateTotalDepth(const std::string& bamFile, const std::string& chr, int position) {
    int totalDepth = 0;
    BamReader newreader(bamFile); 

    std::string region = chr + ":" + std::to_string(position) + "-" + std::to_string(position);
    newreader.SetRegion(region);

    BamRecord record;
    while (newreader.GetNextRecord(record)) {
        totalDepth++;
    }

    return totalDepth;
}

std::string to_uppercase(const std::string& input) {
    std::string result = input;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}

// Function to derive consensus sequences from reads covering a region
std::vector<std::string> deriveConsensusSequences(const std::vector<clustered_aln_t>& reads) {
    std::vector<std::string> contigs;
    std::string currentContig;
    size_t refStart = std::numeric_limits<size_t>::max();
    size_t refEnd = 0;

    // Determine reference span covered by reads
    for (const auto& read : reads) {
        refStart = std::min(refStart, static_cast<size_t>(read.pos));
        refEnd = std::max(refEnd, static_cast<size_t>(read.pos + read.read.Seq().length()));
    }
    // std::cout << "chr1"<< ":" << refStart <<"-" << refEnd << std::endl;
    // Create a pileup for the region
    std::vector<std::unordered_map<char, int>> pileup(refEnd - refStart);

    // Populate the pileup with base counts
    for (const auto& read : reads) {
        const std::string& seq = read.read.Seq();
        size_t start_pos = read.pos - refStart;

        for (size_t i = 0; i < seq.length(); ++i) {
            char base = seq[i];
            pileup[start_pos + i][base]++;
        }
    }

     // Derive consensus sequences
    for (size_t i = 0; i < pileup.size(); ++i) {
        if (!pileup[i].empty()) {
            // Check if depth (number of reads covering this position) is at least 2
            int depth = 0;
            char consensus_base = 'N';  // Default to 'N' (ambiguous) if no votes

            // Calculate depth (total count of reads covering this position)
            for (const auto& [base, count] : pileup[i]) {
                depth += count;
                if (count > pileup[i][consensus_base]) {
                    consensus_base = base;
                }
            }

            // Ensure depth is at least 2 before appending to consensus
            if (depth >= 5) {
                currentContig += consensus_base;  // Append consensus base to current contig
            } else {
                // No coverage at this position, start a new contig if current one is not empty
                if (!currentContig.empty()) {
                    contigs.push_back(currentContig);
                    currentContig.clear();
                }
            }
        } else {
            // No reads covering this position, start a new contig
            if (!currentContig.empty()) {
                contigs.push_back(currentContig);
                currentContig.clear();
            }
        }
    }

    // Finalize the last contig
    if (!currentContig.empty()) {
        contigs.push_back(currentContig);
    }
    return contigs;
}

std::unordered_map<std::string, std::vector<size_t>> buildKmerIndex(const std::string& sequence, size_t k) {
    std::unordered_map<std::string, std::vector<size_t>> kmerIndex;
    size_t n = sequence.length();
    for (size_t i = 0; i <= n - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        kmerIndex[kmer].push_back(i);
    }
    return kmerIndex;
}

std::unordered_map<std::string, int> parseFastaIndex(const std::string& fastaIndexFile) {
    std::unordered_map<std::string, int> chromosomeLengths;
    std::ifstream file(fastaIndexFile);
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string chromosome;
        int length;
        if (iss >> chromosome >> length) {
            chromosomeLengths[chromosome] = length;
        }
    }
    return chromosomeLengths;
}



std::vector<variant_t> DetectComplexIndels(std::map<std::string, std::vector<clustered_aln_t>>& clusteredReads,
    const std::string& bamFile, const std::string& refFasta, const int& threads, const std::string& chromosomeName) {

    BamReader reader(bamFile);
    ThreadSafeQueue<variant_t> uniqueVariants;


    std::string fastaIdx = refFasta + ".fai";
    auto chromosomeLengths = parseFastaIndex(fastaIdx);
    // Set the number of threads for OpenMP
    omp_set_num_threads(threads);

    #pragma omp parallel
    {
        std::vector<variant_t> localVariants; 

        #pragma omp for schedule(dynamic, 1)
        for (size_t idx = 0; idx < clusteredReads.size(); ++idx) {
            auto it = std::next(clusteredReads.begin(), idx);
            const auto& clusterId = *it;

            // Skip clusters not belonging to the current chromosome
            if (clusterId.second.empty() || clusterId.second[0].chromosome != chromosomeName) {
                continue;
            }

            const std::vector<clustered_aln_t>& reads = clusterId.second;
            std::vector<std::string> tmpReads;

            long int refStart = 90000000000;
            long int refEnd = 0;
            std::string chr;

            double meanBaseQuality = 0.00;
            int totalReads = 0;
            for (const auto& read : reads) {
                const std::string& seq = read.read.Seq();
                tmpReads.push_back(seq);
                meanBaseQuality += read.mean_bq;
                totalReads++;
                refStart = std::min(refStart, static_cast<long>(read.pos));
                refEnd = std::max(refEnd, static_cast<long>(read.pos + seq.length()));
                chr = read.chromosome;
            }

            if (totalReads > 0) {
                meanBaseQuality /= totalReads;
            } else {
                meanBaseQuality = 0.0;
            }

            int chromosomeLength = chromosomeLengths[chr];
            RefFasta ref(refFasta);
            refStart = std::max(1L, refStart - 49);  // Adjust for 1-based indexing
            refEnd = std::min(static_cast<long>(chromosomeLength), refEnd + 50);
            std::cout << chr << " " << refStart << " " << refEnd << std::endl;

            std::string regionalFasta = ref.fetchSequence(chr, refStart - 49, refEnd + 50);
            regionalFasta = to_uppercase(regionalFasta);

            const int kmerSize = 20;
            auto kmerIndex = buildKmerIndex(regionalFasta, kmerSize);

            int kSize = 31;
            int maxMismatches = 1;
            bool debug = false;

            Assembler assembler(tmpReads, kSize, maxMismatches, debug);
            vector<string> contigs = assembler.ungappedGreedy();

            std::vector<contig_t> contigRecVector;
            std::vector<int> mapQualVector;

            std::unordered_set<std::string> usedReads;

            for (const auto& contig : contigs) {
                std::string revContig = reverseComplement(contig);

                contig_t contigRecord;
                contigRecord.chr = chr;
                contigRecord.seq = contig;

                long int contigStart = 90000000000;
                long int contigEnd = 0;
                bool isGood = false;
                int readSupport = 0;

                int minusStrand = 0;
                int plusStrand = 0;

                
                for (const auto& read : reads) {

                    if (usedReads.find(read.read.Seq()) != usedReads.end()) {
                        // Skip this read if it has already been used
                        continue;
                    }

                    auto [editDistance, cigar] = alignReadToContig(read.read.Seq(), contig);
                    if (editDistance >= 0 && editDistance < 5) {
                        readSupport++;
                        isGood = true;
                        mapQualVector.push_back(read.read.mapQual());
                        contigStart = std::min(contigStart, static_cast<long>(read.pos));
                        contigEnd = std::max(contigEnd, static_cast<long>(read.pos + read.read.Seq().size()));
                        contigRecord.seq = contig;

                        if (read.strand == '+') {
                            plusStrand++;
                        } else {
                            minusStrand++;
                        }
                         usedReads.insert(read.read.Seq());

                    } 
                    else {
                        std::tie(editDistance, cigar) = alignReadToContig(read.read.Seq(), revContig);
                        if (editDistance >= 0 && editDistance < 5) {
                            readSupport++;
                            isGood = true;
                            mapQualVector.push_back(read.read.mapQual());
                            contigStart = std::min(contigStart, static_cast<long>(read.pos));
                            contigEnd = std::max(contigEnd, static_cast<long>(read.pos + read.read.Seq().size()));
                            contigRecord.seq = revContig;
                            if (read.strand == '+') {
                                plusStrand++;
                            } else {
                                minusStrand++;
                            }
                            usedReads.insert(read.read.Seq());
                        }
                    }
                }
                if (isGood) {
                    contigRecord.start = contigStart;
                    contigRecord.end = contigEnd;
                    contigRecord.readSupport = readSupport;
                    contigRecord.plusStrand = plusStrand;
                    contigRecord.minusStrand = minusStrand;
                    contigRecVector.push_back(contigRecord);
                }
            }

            for (auto& contig : contigRecVector) {
                int match = 2;
                int mismatch = -6;
                int gap_open = -18;
                int gap_extend = 0;
                int match_special = 0;
                int mismatch_special = -20;

                AlignmentResult aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);

                std::string alnCigar = aln.cigarExtended;
                int numOperations = aln.non_match_operations;

                if (numOperations > 3) {
                    match_special = -10;
                    aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
                    alnCigar = aln.cigarExtended;
                    numOperations = aln.non_match_operations;
                }

                int leadingMatches = 0;
                int trailingMatches = 0;
                bool inLeadingSegment = true;

                for (char op : alnCigar) {
                    if (op == 'M') {
                        if (inLeadingSegment) {
                            leadingMatches++;
                        } else {
                            trailingMatches++;
                        }
                    } else {
                        inLeadingSegment = false;
                    }
                }

                if (leadingMatches < 20 || trailingMatches < 20) {
                    continue;
                }

                std::string leadingSeq = aln.seq2_align.substr(0, 20);
                std::string trailingSeq = aln.seq2_align.substr(aln.seq2_align.length() - 20);

                if (kmerIndex[leadingSeq].size() > 1 || kmerIndex[trailingSeq].size() > 1) {
                    continue;
                }

                std::vector<MutationSegment> mutationSegments;
                bool flag = false;
                int start = -1;
                int end = -1;
                const int maxDistance = 10;

                int numOp = 0;
                char prevOp;

                int mod = 0;
                for (size_t i = 0; i < aln.cigarExtended.length(); ++i) {
                    char op = aln.cigarExtended[i];

                    if (!flag && (op == 'X' || op == 'D' || op == 'I')) {
                        flag = true;
                        start = i - mod;
                        end = i - mod;
                        numOp = 1;
                        prevOp = op;
                        mod = 0;
                    } else if (flag) {
                        if (op == 'X' || op == 'D' || op == 'I') {
                            if (i - end <= maxDistance) {
                                end = i;
                                mod = 0;
                            } else {
                                if (numOp > 1 || (numOp == 1 && prevOp != 'X')) {
                                    mutationSegments.push_back({start, end});
                                    start = i;
                                    end = i;
                                }
                                numOp = 1;
                                flag = false;
                                mod = 1;
                            }
                            numOp++;
                        } else {
                            flag = false;
                        }
                    }
                }

                if (numOp > 1 || (numOp == 1 && aln.cigarExtended[start] != 'X')) {
                    mutationSegments.push_back({start, end});
                }

                for (const auto& seg : mutationSegments) {
                    int newstart = seg.start;
                    int newend = seg.end;
                    std::string refAllele;
                    std::string altAllele;

                    if (newstart < 0 || numOperations > 4 || ((newend - newstart) < 2 && numOperations > 3)) {
                        continue;
                    }

                    std::string reducedCigar = alnCigar.substr(newstart - 1, newend + 1 - newstart + 1);
                    int i = aln.query_start + newstart;
                    int j = newstart - 1;

                    for (auto& op : reducedCigar) {
                        char ref_ntd = regionalFasta[i];
                        char alt_ntd = aln.seq2_align[j];
                        if (op != 'D' && op != 'I') {
                            refAllele += ref_ntd;

                            if (alt_ntd != '-') {
                                altAllele += alt_ntd;
                            }
                            i++;
                            j++;
                        } else {
                            if (op == 'D') {
                                refAllele += ref_ntd;
                                i++;
                            } else if (op == 'I') {
                                if (alt_ntd != '-') {
                                    altAllele += alt_ntd;
                                }
                                j++;
                            }
                        }
                    }

                    variant_t variant;
                    long newPos = aln.query_start + refStart + newstart - 49;

                    variant.chr = contig.chr;
                    variant.pos = newPos;
                    variant.ref = refAllele;
                    variant.alt = altAllele;
                    variant.readSupport = contig.readSupport;
                    variant.mapQual = calculateMean(mapQualVector);
                    variant.totalDepth = CalculateTotalDepth(bamFile, variant.chr, variant.pos);
                    int refDepth = variant.totalDepth - variant.readSupport;
                    variant.alleleDepth = std::to_string(refDepth) + "," + std::to_string(variant.readSupport);
                    variant.alleleFrequency = static_cast<float>(variant.readSupport) / variant.totalDepth;
                    variant.kmerDiversity = calculateKmerDiversity(contig.seq, 3);

                    auto [upstreamDiversity, downstreamDiversity] = getFlankingKmerDiversity(variant.chr, variant.pos, refFasta, 50);
                    variant.flankingKmerDiversityUpstream = upstreamDiversity;
                    variant.flankingKmerDiversityDownstream = downstreamDiversity;
                    variant.gcContent = calculateGCContent(contig.seq);

                    auto [longestHomopolymer, numberOfHomopolymers] = calculateHomopolymerMetrics(contig.seq);
                    variant.longestHomopolymerRun = longestHomopolymer;
                    variant.numberOfHomopolymerRuns = numberOfHomopolymers;
                    variant.meanBaseQuality = meanBaseQuality;

                    auto [longestDintd, numberOfDintd] = calculateDinucleotideMetrics(contig.seq);
                    variant.longestDintdRun = longestDintd;
                    variant.numberOfDintd = numberOfDintd;

                    variant.plusStrand = contig.plusStrand;
                    variant.minusStrand = contig.minusStrand;
                    variant.strandBias = ComputeStrandBias(contig.plusStrand, contig.minusStrand);

                    uniqueVariants.push(variant);
                }
            }
        }
    }

    std::vector<variant_t> result;
    variant_t var;
    while (uniqueVariants.try_pop(var)) {
        result.push_back(var);
    }

    return result;
}

