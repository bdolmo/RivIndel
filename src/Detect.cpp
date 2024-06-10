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

double ComputeStrandBias(int forwardStrandCount, int reverseStrandCount) {
    if (forwardStrandCount + reverseStrandCount == 0) {
        return 0.0;
    }
    
    return static_cast<double>(std::fabs(forwardStrandCount - reverseStrandCount)) / (forwardStrandCount + reverseStrandCount);
}


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

void PopulateVCF(const std::map<std::string, variant_t>& variants, const std::string& bamFile, const std::string& vcfFilename) {

    std::vector<variant_t> allVariants;
    for (const auto& kv : variants) {
        allVariants.push_back(kv.second);
    }

    std::sort(allVariants.begin(), allVariants.end(), CompareVariants);

    VcfWriter writer(vcfFilename, bamFile);
    writer.writeHeader();
    for (const auto& var : allVariants) {
        writer.writeVariant(var);
    }
}


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


std::map<std::string, variant_t> DetectComplexIndels(std::map<std::string, std::vector<clustered_aln_t>>& clusteredReads,
    const std::string& bamFile, const std::string& refFasta, const int& threads) {
    RefFasta ref(refFasta);
    BamReader reader(bamFile);
    std::map<std::string, variant_t> uniqueVariants;

    // Set the number of threads for OpenMP
    omp_set_num_threads(threads);

    #pragma omp parallel
    {

        #pragma omp for schedule(dynamic, 1)
        for (size_t idx = 0; idx < clusteredReads.size(); ++idx) {
            auto it = std::next(clusteredReads.begin(), idx);
            const auto& clusterId = *it;

            std::string msg = " INFO: Analyzing cluster " + clusterId.first;
            std::cout << msg << std::endl;

            const std::vector<clustered_aln_t>& reads = clusterId.second;
            std::vector<std::string> tmpReads;

            long int refStart = 90000000000;
            long int refEnd = 0;
            std::string chr;

            double meanBaseQuality = 0.00;
            std::vector<double> mbqVec;
            int totalReads = 0;
            for (const auto& read : reads) {
                const std::string& seq = read.read.Seq();
                tmpReads.push_back(seq);
                // std::cout << read.mean_bq << std::endl;
                meanBaseQuality += read.mean_bq; // Sum up the mean base quality values
                totalReads++;
                refStart = std::min(refStart, static_cast<long>(read.pos));
                refEnd = std::max(refEnd, static_cast<long>(read.pos + seq.length()));
                chr = read.chromosome;
                // std::cout << clusterId.first << " " << read.chromosome << ":" << read.pos << " " << read.strand<< std::endl;
            }

            if (totalReads > 0) {
                meanBaseQuality /= totalReads; // Calculate the average mean base quality
            } 
            else {
                meanBaseQuality = 0.0; // Handle the case where there are no reads
            }

            std::string regionalFasta = ref.fetchSequence(chr, refStart - 49, refEnd + 50);
            regionalFasta = to_uppercase(regionalFasta);
            std::vector<std::string> contigs = assembleWithFermiLite(tmpReads);

            std::vector<contig_t> contigRecVector;
            std::vector<int> mapQualVector;

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
                        }
                        else {
                            minusStrand++;
                        }
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
                            }
                            else {
                                minusStrand++;
                            }
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

                // std::cout << "before: " << contig.start << std::endl;
                // contig.start = aln.ref_start+refStart;
                // std::cout << "refStart:" << refStart << " " << aln.query_start << " " << contig.start << std::endl;

                std::string alnCigar = aln.cigarExtended;
                int numOperations = aln.non_match_operations;

                if (numOperations > 3) {
                    match_special = -10;
                    aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
                    alnCigar = aln.cigarExtended;
                    numOperations = aln.non_match_operations;
                }

                // Check for at least 20 flanking matches
                int leadingMatches = 0;
                int trailingMatches = 0;
                bool inLeadingSegment = true;
                std::string leadingSeq;
                std::string trailingSeq;

                for (char op : alnCigar) {
                    if (op == 'M') {
                        if (inLeadingSegment) {
                            leadingMatches++;
                            if (leadingMatches <= 20) {
                                leadingSeq += contig.seq[leadingMatches - 1];
                            }
                        } else {
                            trailingMatches++;
                            if (trailingMatches <= 20) {
                                trailingSeq += contig.seq[contig.seq.size() - trailingMatches];
                            }
                        }
                    } else {
                        inLeadingSegment = false; 
                    }
                }
                if (leadingMatches < 20 || trailingMatches < 20) {
                    continue;
                }

                std::vector<std::pair<int, int>> mutationSegments;
                bool flag = false;
                int start = -1;
                int end = -1;
                const int maxDistance = 10;

                int numOp = 0;
                char prevOp;

                int mod = 0;
                for (size_t i = 0; i < alnCigar.length(); ++i) {
                    char op = alnCigar[i];
                                    
                // MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM X MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDDDDDDDMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
                // ATAAGAAGTGTCATAGGTTCGCGTGTTTTTCACAAC T GATTATGATGTCCATGATGACAGCGTGCACTAAGCGCTGTGTCCAACCAACCCTTTACCTGGCATTTTAGCCAGCAAGACCCGCTTGCATTGCAGTCGTTCACGGTGCTCGCTGGCACAATTTTAGAATGTTACATGAA
                // |||||||||||||||||||||||||||||||||||| | ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||       ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                // ATAAGAAGTGTCATAGGTTCGCGTGTTTTTCACAAC a G ATTATGATGTCCATGATGACAGCGTGCACTAAGCGCTGTGTCCAACCAACCCTTTACCT-------TAGCCAGCAAGACCCGCTTGCATTGCAGTCGTTCACGGTGCTCGCTGGCACAATTTTAGAATGTTACATGAA

                    if (!flag && (op == 'X' || op == 'D' || op == 'I')) {
                        flag = true;
                        start = i-mod;
                        end = i-mod;
                        numOp = 1;
                        prevOp = op;
                        mod = 0;
                    } 
                    else if (flag) {
                        if (op == 'X' || op == 'D' || op == 'I') {
                            if (i - end <= maxDistance) {
                                end = i;
                                mod = 0;
                            }
                            else {
                                if (numOp > 1 || (numOp == 1 && prevOp != 'X')) {
                                    mutationSegments.push_back({start, end});
                                    start = i;
                                    end = i;
                                }
                                numOp = 1;
                                flag  = false;
                                mod = 1;

                            }
                            numOp++;
                        }
                        else {
                            flag = false;
                        }
                    }
                }

                if (numOp > 1 || (numOp == 1 && alnCigar[start] != 'X')) {
                    mutationSegments.push_back({start, end});
                }


                std::cout << start << " " << end << " " << numOperations << std::endl;

                for (const auto& seg : mutationSegments) {
                    start = seg.first;
                    end = seg.second;

                    if (start >= 0 && numOperations < 5 && (end - start >= 2 || numOperations <= 3)) {
                        std::string refAllele;
                        std::string altAllele;
                        std::string reducedCigar = alnCigar.substr(start - 1, end + 1 - start + 1);
                        int i = aln.query_start + start;
                        int j = start - 1;


                        std::cout << alnCigar << std::endl;
                        std::cout << aln.seq1_align << std::endl;
                        std::cout << aln.spacer << std::endl;
                        std::cout << aln.seq2_align << std::endl;
                        std::cout << aln.query_start << "=>" << aln.query_end << " COORD: " << aln.query_start + refStart <<  std::endl;
                        std::cout << aln.ref_start << "=>" << aln.ref_end << std::endl;
                        std::cout << reducedCigar << std::endl;
                        
                        for (auto& op : reducedCigar) {
                            // char ref_ntd = aln.seq1_align[i];
                            char ref_ntd = regionalFasta[i];
                            char alt_ntd = aln.seq2_align[j];
                            if (op != 'D' && op != 'I') {
                                refAllele += ref_ntd;

                                if (alt_ntd != '-') {
                                    altAllele += alt_ntd;
                                }
                                i++;
                                j++;
                            }
                            else {
                                if (op == 'D') {
                                    refAllele += ref_ntd;
                                    i++;
                                }
                                else if (op == 'I') {
                                    if (alt_ntd != '-') {
                                        altAllele += alt_ntd;
                                    }
                                    // if (ref_ntd != '-') {
                                    //     refAllele += ref_ntd;
                                    // }
                                    j++;
                                }
                                // if (op == 'X') {
                                //     j++;
                                // }
                            }

                        }
                        std::string variantKey = contig.chr + ":" + std::to_string(aln.query_start + refStart+start-49) + ":" + refAllele + ">" + altAllele;
                        #pragma omp critical
                        {
                            if (uniqueVariants.find(variantKey) != uniqueVariants.end()) {
                                uniqueVariants[variantKey].readSupport++;
                            } 
                            else {
                                std::cout << "CONTIG: "<< contig.chr << " " << contig.start  <<  " " << contig.end << std::endl;
                                variant_t variant;

                                long newPos = aln.query_start + refStart+start-49;

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
                                
                                // Calculate flanking k-mer diversity
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

                                uniqueVariants[variantKey] = variant;
                                std::cout << variant.chr << ":" << variant.pos << " " << variant.ref << ">" << variant.alt << std::endl << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    return uniqueVariants;
}

