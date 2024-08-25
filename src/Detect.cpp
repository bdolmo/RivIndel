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
#include <regex>
#include "Assembler.h"
#include <unordered_set>
#include "ssw_cpp.h"


AlignmentResult align_with_ssw(const std::string& read, const std::string& contig) {
    // Initialize a default scoring matrix for DNA sequences
    const int8_t match = 2;     // Match score
    const int8_t mismatch = 8;  // Mismatch penalty
    const int8_t gapOpen = 10;  // Gap open penalty
    const int8_t gapExtend = 1; // Gap extend penalty

    // Initialize the SSW aligner
    StripedSmithWaterman::Aligner aligner(match, mismatch, gapOpen, gapExtend);

    // Initialize the SSW filter (optional settings, can be adjusted as needed)
    StripedSmithWaterman::Filter filter;

    // Create an alignment object to store the result
    StripedSmithWaterman::Alignment alignment;

    // Perform the alignment
    aligner.Align(read.c_str(), contig.c_str(), contig.size(), filter, &alignment);

    // Extract the compacted CIGAR string
    std::string compactedCigar = alignment.cigar_string;

    // Generate the extended CIGAR string and adjust query_start and query_end to skip soft-clips
    std::string cigarExtended;
    int query_start = alignment.query_begin;
    int query_end = alignment.query_end;
    int ref_start = alignment.ref_begin;
    int ref_end = alignment.ref_end;

    // Count the number of non-'M' operations from the compacted CIGAR string
    std::set<char> operations;
    for (size_t i = 0; i < compactedCigar.length(); ++i) {
        char c = compactedCigar[i];
        if (std::isalpha(c) && c != 'M' && c != 'S') {
            operations.insert(c);  // Insert into set to ensure uniqueness
        }
    }

    // Count of distinct non-'M' and non-'S' operations
    int non_match_operations = operations.size();

    std::string seq1_align, seq2_align, spacer;
    // Optionally, generate the extended CIGAR if needed
    int readPos = alignment.query_begin;
    int contigPos = alignment.ref_begin;
    for (size_t i = 0; i < compactedCigar.length(); ++i) {
        char op = compactedCigar[i];
        if (isdigit(op)) {
            int len = 0;
            while (isdigit(compactedCigar[i])) {
                len = len * 10 + (compactedCigar[i] - '0');
                ++i;
            }
            op = compactedCigar[i];
            if (op == 'S') {
                continue;
            }

            for (int j = 0; j < len; ++j) {
                cigarExtended += op;
                if (op == 'M') { // Match/mismatch
                    seq1_align += read[readPos];
                    seq2_align += contig[contigPos];
                    spacer += (read[readPos] == contig[contigPos] ? '|' : '*'); // Match or mismatch
                    readPos++;
                    contigPos++;
                } else if (op == 'I') { // Insertion in the read
                    seq1_align += read[readPos];
                    seq2_align += '-';
                    spacer += ' ';
                    readPos++;
                    // non_match_operations++;  // Count as non-match operation
                } else if (op == 'D') { // Deletion in the read (gap in the read)
                    seq1_align += '-';
                    seq2_align += contig[contigPos];
                    spacer += ' ';
                    contigPos++;
                    // non_match_operations++;  // Count as non-match operation
                }
                // } else if (op == 'S') { // Soft clipping
                //     seq1_align += read[readPos];  // Include soft-clipped bases in the read alignment
                //     seq2_align += '-';            // No reference base, so align to a gap
                //     spacer += ' ';                // No match or mismatch indication
                //     readPos++;  // Move forward in the read but not in the reference
                //     // contigPos++;
                // }
            }
            // for (int j = 0; j < len; ++j) {
            //     cigarExtended += op;
            //     if (op == 'M') { // Match/mismatch
            //         seq1_align += read[readPos];
            //         seq2_align += contig[contigPos];
            //         spacer += (read[readPos] == contig[contigPos] ? '|' : '*'); // Match or mismatch
            //         readPos++;
            //         contigPos++;
            //     } else if (op == 'I') { // Insertion in the read
            //         seq1_align += read[readPos];
            //         seq2_align += '-';
            //         spacer += ' ';
            //         readPos++;
            //         non_match_operations++;  // Count as non-match operation
            //     } else if (op == 'D') { // Deletion in the read (gap in the read)
            //         seq1_align += '-';
            //         seq2_align += contig[contigPos];
            //         spacer += ' ';
            //         contigPos++;
            //         non_match_operations++;  // Count as non-match operation
            //     }
            // }
        }
    }

    // std::cout << seq1_align << std::endl;
    // std::cout << spacer << std::endl;
    // std::cout << seq2_align << std::endl<< std::endl;

    // Create the AlignmentResult2 struct
    AlignmentResult result;
    result.query_start = ref_start-1;
    result.query_end = ref_end+1;
    result.ref_start = query_start;
    result.ref_end = query_end+1;
    result.compactedCigar = compactedCigar;
    result.cigarExtended = cigarExtended;
    result.seq1_align = seq1_align;
    result.seq2_align = seq2_align;
    result.spacer = spacer;

    result.non_match_operations = non_match_operations;  // Store the count of non-'M' operations

    return result;
}

// Function to check for clusters of different operations within 10 positions
bool hasOpCluster(const std::string& cigar, const int& numOp) {
    // Sliding window to store operations
    for (size_t i = 0; i <= cigar.size() - 10; ++i) {
        std::unordered_set<char> uniqueOps;
        for (size_t j = i; j < i + 10; ++j) {
            uniqueOps.insert(cigar[j]);
        }
        if (uniqueOps.size() >= 3) {
            return true;
        }
    }
    return false;
}

// Function to extract k-mers from a sequence
std::vector<std::string> extractKmers(const std::string& sequence, int k) {
    std::vector<std::string> kmers;
    if (sequence.size() < k) {
        return kmers;  // Return an empty vector if sequence is too short
    }
    kmers.reserve(sequence.size() - k + 1);
    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        kmers.emplace_back(sequence.substr(i, k));
    }
    return kmers;
}


// Function to count k-mer frequencies
std::unordered_map<std::string, int> countKmerFrequencies(const std::vector<std::string>& sequences, int k) {
    std::unordered_map<std::string, int> kmerCounts;
    for (const auto& sequence : sequences) {
        if (sequence.size() < k) {
            continue;  // Skip sequences that are too short
        }
        std::vector<std::string> kmers = extractKmers(sequence, k);
        for (const auto& kmer : kmers) {
            kmerCounts[kmer]++;
        }
    }
    return kmerCounts;
}


// Function to identify and correct erroneous k-mers
std::string correctSequence(const std::string& sequence, int k, const std::unordered_map<std::string, int>& kmerCounts, int threshold) {
    if (sequence.length() < k) {
        return sequence;  // Return the sequence as is if it's too short to have k-mers
    }
    std::string correctedSequence = sequence;
    for (size_t i = 0; i <= sequence.length() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        if (kmerCounts.at(kmer) < threshold) {
            std::string bestKmer = kmer;
            int bestCount = 0;
            for (size_t j = 0; j < k; ++j) {
                std::string tempKmer = kmer;
                for (char nucleotide : {'A', 'C', 'G', 'T'}) {
                    tempKmer[j] = nucleotide;
                    int tempCount = kmerCounts.count(tempKmer) ? kmerCounts.at(tempKmer) : 0;
                    if (tempCount > bestCount) {
                        bestCount = tempCount;
                        bestKmer = tempKmer;
                    }
                }
            }
            for (size_t j = 0; j < k; ++j) {
                correctedSequence[i + j] = bestKmer[j];
            }
        }
    }
    return correctedSequence;
}


// Main function to correct a vector of sequences
std::vector<std::string> correctSequences(const std::vector<std::string>& sequences, int k, int threshold) {
    std::vector<std::string> correctedSequences;
    correctedSequences.reserve(sequences.size());
    std::unordered_map<std::string, int> kmerCounts = countKmerFrequencies(sequences, k);
    
    for (const auto& sequence : sequences) {
        if (sequence.size() < k) {
            correctedSequences.push_back(sequence);  // Add the sequence unchanged if it's too short
        } else {
            correctedSequences.push_back(correctSequence(sequence, k, kmerCounts, threshold));
        }
    }
    
    return correctedSequences;
}


double ComputeStrandBias(int forwardStrandCount, int reverseStrandCount) {
    if (forwardStrandCount + reverseStrandCount == 0) {
        return 0.0;
    }
    
    return static_cast<double>(std::fabs(forwardStrandCount - reverseStrandCount)) / (forwardStrandCount + reverseStrandCount);
}

struct contig_t{
    std::string seq;
    std::string chr;
    int start;
    int end;
    int readSupport;
    int plusStrand;
    int minusStrand;
    std::vector<clustered_aln_t> contigReads;
};

struct Locus {
    std::string chromosome;
    int start;
    int end;
    std::vector<contig_t> contigs;
};

// struct Contig {
//     std::string chromosome;
//     int start;
//     int end;
//     int pos;
//     std::string seq;
//     std::string cigarExtended;
// }

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
    // Initialize a default scoring matrix for DNA sequences
    const int8_t match = 2;      // Match score
    const int8_t mismatch = 8;   // Mismatch penalty
    const int8_t gapOpen = 10;   // Gap open penalty
    const int8_t gapExtend = 1;  // Gap extend penalty

    // Initialize the SSW aligner
    StripedSmithWaterman::Aligner aligner(match, mismatch, gapOpen, gapExtend);

    // Initialize the SSW filter (optional settings, can be adjusted as needed)
    StripedSmithWaterman::Filter filter;

    // Create an alignment object to store the result
    StripedSmithWaterman::Alignment alignment;

    // Perform the alignment
    aligner.Align(read.c_str(), contig.c_str(), contig.size(), filter, &alignment);

    // Calculate the aligned length
    int alignedLength = alignment.ref_end - alignment.ref_begin + 1;

    // Check if the aligned length is less than 15 bases
    if (alignedLength < 15) {
        return {-1, ""};  // Return -1 for mismatches and an empty CIGAR if not valid
    }

    // Check for valid 5' or 3' end alignment
    bool alignsAt5Prime = (alignment.query_begin == 0 || alignment.ref_begin == 0);
    bool alignsAt3Prime = (alignment.query_end == read.size() - 1 || alignment.ref_end == contig.size() - 1);

    // If the alignment does not cover at least one end, it's not valid
    if (!alignsAt5Prime && !alignsAt3Prime) {
        return {-1, ""};  // Return -1 if it doesn't align to either end
    }

    // Convert the CIGAR vector to a string
    std::string cigar = alignment.cigar_string;

    // Check if the alignment is continuous (i.e., no large gaps or soft clips in the overlap regions)
    bool isContinuous = true;
    int currentPos = 0;
    bool hasSoftClip = false;
    for (size_t i = 0; i < cigar.size(); ++i) {
        char op = cigar[i];
        if (isdigit(op)) {
            // Extract the length of the operation
            int len = 0;
            while (isdigit(cigar[i])) {
                len = len * 10 + (cigar[i] - '0');
                ++i;
            }
            op = cigar[i];  // Get the operation character (M, I, D, S, etc.)

            // Check for insertions, deletions, and soft-clips
            if (op == 'I' || op == 'D' || op == 'S') {
                // If the operation occurs near the start or end of the alignment and is significant, discard it
                if ((currentPos <= 10 || currentPos >= cigar.size() - 10) && len >= 5) {
                    isContinuous = false;
                    break;
                }
                // Check for soft-clipped regions, discard if found
                if (op == 'S' && len >= 10) {
                    hasSoftClip = true;
                    break;  // Soft-clipped regions are not allowed in overlap alignments
                }
            }
            currentPos += len;
        }
    }

    // If the alignment has large gaps or soft-clips in the overlap, discard it
    if (!isContinuous || hasSoftClip) {
        return {-1, ""};  // Discard alignments with soft-clips or non-continuous overlaps
    }

    // Calculate the total number of mismatches
    int totalMismatches = alignment.mismatches;

    // Generate the aligned sequence output
    std::string alignedRead = read.substr(alignment.query_begin, alignedLength);
    std::string alignedContig = contig.substr(alignment.ref_begin, alignedLength);

    // Return the total number of mismatches and the CIGAR string
    return {totalMismatches, cigar};
}
// std::pair<int, std::string> alignReadToContig(const std::string& read, const std::string& contig) {
//     EdlibAlignResult result = edlibAlign(read.c_str(), read.size(), contig.c_str(), contig.size(),
//         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

//     int editDistance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
//     std::string cigar;

//     if (result.status == EDLIB_STATUS_OK && result.alignment) {
//         char* cigarC = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
//         cigar = std::string(cigarC);
//         free(cigarC);
//     }


//     // Print out the offsets and the CIGAR string
//     if (result.status == EDLIB_STATUS_OK) {
//         // std::cout << "Read Start: " << result.startLocations[0] << std::endl;
//         // std::cout << "Read End: " << result.endLocations[0] << std::endl;
//         // std::cout << "Contig Start: " << result.startLocations[1] << std::endl;
//         // std::cout << "Contig End: " << result.endLocations[1] << std::endl;
//         // std::cout << read << " " << contig << std::endl;
//     // printf("%d\n", result.editDistance);
//     // printf("%d\n", result.alignmentLength);
//     // printf("%d\n", result.endLocations[0]);
//         // std::cout << result.alignmentLength << std::endl;

//         // std::cout << "Editdistance: " << editDistance << " CIGAR String: " << cigar << " " << result.alignment << std::endl << std::endl;
//     }

//     edlibFreeAlignResult(result);
//     return {editDistance, cigar};
// }


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

// Helper function to calculate median
double CalculateMedian(std::vector<int>& depths) {
    if (depths.empty()) return 0.0;
    std::sort(depths.begin(), depths.end());
    size_t size = depths.size();
    if (size % 2 == 0) {
        return (depths[size / 2 - 1] + depths[size / 2]) / 2.0;
    } else {
        return depths[size / 2];
    }
}

struct AlignmentStats {
    double errorRate;
    double chimericRate;
    double softClippedRate;
    double medianCoverage;
};

// Function to calculate various alignment statistics including median coverage
AlignmentStats CalculateAlignmentStats(const std::string& bamFile, const std::string& chr, int position) {
    int totalBases = 0;
    int mismatches = 0;
    int totalReads = 0;
    int chimericReads = 0;
    int softClippedReads = 0;

    BamReader newreader(bamFile);

    // Set the region for the desired window
    std::string region = chr + ":" + std::to_string(position - 100) + "-" + std::to_string(position + 100);
    newreader.SetRegion(region);

    int maxReads = 500;
    int count = 0;
    BamRecord record;

    // Use a map to track the depth at each position
    std::map<int, int> depthMap;

    while (newreader.GetNextRecord(record)) {
        totalReads++;
        std::string mdTag = record.MDtag();
        totalBases += record.Seq().length();
        count++;

        // Track the read depth at the positions spanned by this read
        long readStart = record.Position();
        long readEnd = readStart + record.Seq().length();
        for (long pos = readStart; pos < readEnd; ++pos) {
            if (pos >= position - 10 && pos <= position) {  // Only count positions in the window
                depthMap[pos]++;
            }
        }

        // Check for chimeric reads (mapped to different chromosomes or far apart)
        if (record.chrName() != record.chrMateName()) {
            chimericReads++;
        }

        // Check for soft-clipped reads
        std::string cigar = record.cigarString();
        bool hasSoftClip = cigar.find('S') != std::string::npos;
        if (hasSoftClip) {
            softClippedReads++;
        }

        if (count > maxReads) {
            break;
        }

        bool foundDel = false;
        // Manually parse the MD tag to find mismatches
        for (size_t i = 0; i < mdTag.size(); ++i) {
            if (std::isdigit(mdTag[i])) {
                if (foundDel) {
                    foundDel = false;
                }

                // Skip the number of matched bases
                while (i < mdTag.size() && std::isdigit(mdTag[i])) {
                    ++i;
                }
                --i;  // Correct the position after exiting the while loop
            } else {
                // Count mismatches (non-numeric parts of the MD tag)
                if (mdTag[i] == '^') {
                    foundDel = true;
                    continue;
                }
                if (!foundDel && std::isalpha(mdTag[i])) {
                    mismatches++;
                }
            }
        }
    }

    // Collect depths from the map for the range [position-10, position]
    std::vector<int> depths;
    for (int pos = position - 10; pos <= position; ++pos) {
        if (depthMap.find(pos) != depthMap.end()) {
            depths.push_back(depthMap[pos]);
        } else {
            depths.push_back(0);  // If no reads cover this position, assume depth of 0
        }
    }

    // Calculate the median depth
    double medianDepth = CalculateMedian(depths);

    AlignmentStats stats;
    if (totalBases == 0) {
        stats.errorRate = 0.0;
    } else {
        stats.errorRate = static_cast<double>(mismatches) / totalBases;
    }
    stats.chimericRate = static_cast<double>(chimericReads) / totalReads;
    stats.softClippedRate = static_cast<double>(softClippedReads) / totalReads;
    stats.medianCoverage = medianDepth;

    return stats;
}



// struct AlignmentStats {
//     double errorRate;
//     double chimericRate;
//     double softClippedRate;
// };

// // Function to calculate various alignment statistics
// AlignmentStats CalculateAlignmentStats(const std::string& bamFile, const std::string& chr, int position) {
//     int totalBases = 0;
//     int mismatches = 0;
//     int totalReads = 0;
//     int chimericReads = 0;
//     int softClippedReads = 0;
//     BamReader newreader(bamFile); 

//     std::string region = chr + ":" + std::to_string(position - 100) + "-" + std::to_string(position + 100);
//     newreader.SetRegion(region);
//     int maxReads = 500;
//     int count = 0;
//     BamRecord record;
//     while (newreader.GetNextRecord(record)) {
//         totalReads++;
//         std::string mdTag = record.MDtag();
//         totalBases += record.Seq().length();
//         count++;

//         // Check for chimeric reads (mapped to different chromosomes or far apart)
//         if (record.chrName() != record.chrMateName()) {
//             chimericReads++;
//         }

//         // Check for soft-clipped reads
//         std::string cigar = record.cigarString();
//         bool hasSoftClip = cigar.find('S') != std::string::npos;
//         if (hasSoftClip) {
//             softClippedReads++;
//         }

//         if (count > maxReads) {
//             break;
//         }

//         bool foundDel = false;
//         // Manually parse the MD tag to find mismatches
//         for (size_t i = 0; i < mdTag.size(); ++i) {
//             if (std::isdigit(mdTag[i])) {
//                 if (foundDel) {
//                     foundDel = false;
//                 }

//                 // Skip the number of matched bases
//                 while (i < mdTag.size() && std::isdigit(mdTag[i])) {
//                     ++i;
//                 }
//                 --i; // Correct the position after exiting the while loop
//             } else {
//                 // Count mismatches (non-numeric parts of the MD tag)
//                 if (mdTag[i] == '^') {
//                     foundDel = true;
//                     continue;
//                 }
//                 if (!foundDel && std::isalpha(mdTag[i])) {
//                     mismatches++;
//                 }
//             }
//         }
//     }

//     AlignmentStats stats;
//     if (totalBases == 0) {
//         stats.errorRate = 0.0;
//     } else {
//         stats.errorRate = static_cast<double>(mismatches) / totalBases;
//     }
//     stats.chimericRate = static_cast<double>(chimericReads) / totalReads;

//     stats.softClippedRate = static_cast<double>(softClippedReads) / totalReads;

//     return stats;
// }


// // Function to calculate median
// double CalculateMedian(std::vector<int>& depths) {
//     if (depths.empty()) return 0.0;
//     std::sort(depths.begin(), depths.end());
//     size_t size = depths.size();
//     if (size % 2 == 0) {
//         return (depths[size / 2 - 1] + depths[size / 2]) / 2.0;
//     } else {
//         return depths[size / 2];
//     }
// }

// int CalculateTotalDepthAtPosition(const std::string& bamFile, const std::string& chr, int position) {
//     int depth = 0;
//     BamReader newreader(bamFile); 

//     // Set region for a single position
//     std::string region = chr + ":" + std::to_string(position+1) + "-" + std::to_string(position+1);
//     newreader.SetRegion(region);

//     BamRecord record;
//     while (newreader.GetNextRecord(record)) {
//         depth++;
//     }

//     return depth;
// }

// double CalculateMedianCoverage(const std::string& bamFile, const std::string& chr, int position) {
//     std::vector<int> depths;

//     // Iterate over the window from position-9 to position
//     for (int pos = position - 10; pos <= position-1; pos++) {
//         int depthAtPos = CalculateTotalDepthAtPosition(bamFile, chr, pos);
//         depths.push_back(depthAtPos);  // Store depth at each position
//     }

//     // Calculate and return the median of the depths
//     return CalculateMedian(depths);
// }





int CalculateTotalDepth(const std::string& bamFile, const std::string& chr, int position) {
    int totalDepth = 0;
    BamReader newreader(bamFile); 

    std::string region = chr + ":" + std::to_string(position+1) + "-" + std::to_string(position+1);
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
            } 
            else {
                // No coverage at this position, start a new contig if current one is not empty
                if (!currentContig.empty()) {
                    contigs.push_back(currentContig);
                    currentContig.clear();
                }
            }
        } 
        else {
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


std::vector<std::string> correctReads(std::vector<std::string>& reads) {
    int n_seqs = reads.size();
    bseq1_t *seqs = new bseq1_t[n_seqs]();

    // Populate the bseq1_t array with the input reads
    for (int i = 0; i < n_seqs; ++i) {
        seqs[i].l_seq = reads[i].length();
        seqs[i].seq = strdup(reads[i].c_str());

        // Allocate and set quality scores (all set to 'I' which is ASCII 73)
        seqs[i].qual = new char[reads[i].length() + 1];
        memset(seqs[i].qual, 'I', reads[i].length());
        seqs[i].qual[reads[i].length()] = '\0';
    }

    // Initialize default options
    fml_opt_t opt;
    fml_opt_init(&opt);

    // Perform error correction
    float kcov = fml_correct(&opt, n_seqs, seqs);

    // Optional: adjust parameters based on corrected sequences
    fml_opt_adjust(&opt, n_seqs, seqs);

    // Collect corrected reads into a new vector
    std::vector<std::string> correctedReads;
    for (int i = 0; i < n_seqs; ++i) {
        correctedReads.push_back(std::string(seqs[i].seq));
    }

    // Clean up
    for (int i = 0; i < n_seqs; ++i) {
        free(seqs[i].seq);
        delete[] seqs[i].qual;
    }
    delete[] seqs;

    return correctedReads;
}


bool has_valid_end_operations(const std::string& cigarExtended) {
    int startMatches = 0, startSoftClips = 0;
    int endMatches = 0, endSoftClips = 0;

    // Check the first 5 operations
    for (size_t i = 0; i < cigarExtended.length() && i < 5; ++i) {
        if (cigarExtended[i] == 'M') {
            startMatches++;
        } else if (cigarExtended[i] == 'S') {
            startSoftClips++;
        }
    }

    // Check the last 5 operations
    for (size_t i = cigarExtended.length(); i-- > cigarExtended.length() - 5 && i >= 0; ) {
        if (cigarExtended[i] == 'M') {
            endMatches++;
        } else if (cigarExtended[i] == 'S') {
            endSoftClips++;
        }
    }

    // Valid if there are at least 5 matches or soft clips at both ends
    return (startMatches >= 5 || startSoftClips >= 5) && (endMatches >= 5 || endSoftClips >= 5);
}

std::vector<variant_t> DetectComplexIndels(std::map<std::string, std::vector<clustered_aln_t>>& clusteredReads,
    const std::string& bamFile, const std::string& refFasta, const int& threads, const std::string& chromosomeName, 
    const bool contigsOut, const std::string& contigsOutFile) {

    BamReader reader(bamFile);
    ThreadSafeQueue<variant_t> uniqueVariants;

    ThreadSafeQueue<MutationSegment> mutationSegmentsQueue;

    std::string fastaIdx = refFasta + ".fai";
    auto chromosomeLengths = parseFastaIndex(fastaIdx);
    // Set the number of threads for OpenMP
    omp_set_num_threads(threads);

    #pragma omp parallel
    {
        std::vector<variant_t> localVariants;
        std::vector<Locus> locusData;

        std::unordered_set<std::string> threadLocalUsedReads;
        std::unordered_set<std::string> threadLocalUsedReads2;

        #pragma omp for schedule(dynamic, 1)
        for (size_t idx = 0; idx < clusteredReads.size(); ++idx) {
            auto it = std::next(clusteredReads.begin(), idx);
            const auto& clusterId = *it;

            // Skip clusters not belonging to the current chromosome
            if (clusterId.second.empty() || clusterId.second[0].chromosome != chromosomeName) {
                continue;
            }

            const std::vector<clustered_aln_t>& reads = clusterId.second;
            // std::cout << " TOTAL_READS " << reads.size() << std::endl;

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

            int flankSpace = 49;
            int chromosomeLength = chromosomeLengths[chr];
            RefFasta ref(refFasta);
            refStart = std::max(1L, refStart - flankSpace);  // Adjust for 1-based indexing
            refEnd = std::min(static_cast<long>(chromosomeLength), refEnd + flankSpace+1);
            // std::cout << " INFO: Interrogating "<< chr << ":" << refStart << "-" << refEnd << std::endl;

            Locus newLocus;
            newLocus.chromosome = chr;
            newLocus.start = refStart;
            newLocus.end = refEnd;

            std::string regionalFasta = ref.fetchSequence(chr, refStart - flankSpace, refEnd + flankSpace+1);
            regionalFasta = to_uppercase(regionalFasta);

            const int kmerSize = 20;
            auto kmerIndex = buildKmerIndex(regionalFasta, kmerSize);

            int kSize = 31;
            int maxMismatches = 1;
            bool debug = false;
            std::vector<std::string> corrReads = correctSequences(tmpReads, 21, 2);

            // std::cout << " INPUT CONTIGS" << std::endl;
            Assembler assembler(corrReads, kSize, maxMismatches, debug);

            // Assembler assembler(tmpReads, kSize, maxMismatches, debug);
            vector<string> contigs = assembler.ungappedGreedy();

            // std::cout << "OUTPUT CONTIGS" << std::endl;

            std::vector<contig_t> contigRecVector;
            std::vector<int> mapQualVector;
            std::unordered_set<std::string> usedReads;
            std::sort(contigs.begin(), contigs.end(), [](const std::string& a, const std::string& b) {
                return a.size() > b.size(); // Sort by length, longest strings first
            });

            for (const auto& contig : contigs) {

                // std::cout << "Contig to be aligned with reds: "<< chr << ":" << refStart << "-" << refEnd << " " << contig << std::endl;
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
                std::mutex mtx;
                for (int i=0; i<reads.size();i++) {
                    std::string readName = reads[i].readName;
                    std::string readSeq = reads[i].read.Seq();

                    std::string tag = readSeq + " " + readName;
                    // Use thread-local copy of usedReads to avoid race conditions
                    if (threadLocalUsedReads.find(tag) != threadLocalUsedReads.end()) {
                        continue;  // Skip if already used in this thread
                    }

                    // {
                    //     // Lock the mutex to safely access `usedReads`
                    //     std::lock_guard<std::mutex> lock(mtx);

                    //     // Check if the read has already been used
                    //     if (usedReads.find(tag) != usedReads.end()) {
                    //         continue;  // Skip if already used
                    //     }
                    // }
                    mapQualVector.push_back(reads[i].read.mapQual());
                    auto [editDistance, cigar] = alignReadToContig(readSeq, contig);
                    // std::cout << editDistance << " " << cigar << std::endl<< std::endl;

                    if (editDistance <= 10 && editDistance >= 0) {
                        readSupport++;
                        isGood = true;
                        contigStart = std::min(contigStart, static_cast<long>(reads[i].pos));
                        contigEnd = std::max(contigEnd, static_cast<long>(reads[i].pos + readSeq.size()));
                        contigRecord.seq = contig;
                        contigRecord.contigReads.push_back(reads[i]);

                        if (reads[i].strand == '+') {
                           plusStrand++;
                        } 
                        else {
                           minusStrand++;
                        }
                        threadLocalUsedReads.insert(tag);


                        // {
                        //     // Lock the mutex again to modify `usedReads`
                        //     std::lock_guard<std::mutex> lock(mtx);
                        //     usedReads.insert(tag);
                        // }
                    }
                }
                // std::cout << readSupport << std::endl;
                if (isGood) {
                    contigRecord.start = contigStart;
                    contigRecord.end = contigEnd;
                    contigRecord.readSupport = readSupport;
                    contigRecord.plusStrand = plusStrand;
                    contigRecord.minusStrand = minusStrand;
                    contigRecVector.push_back(contigRecord);
                }
            }

            for (auto& contig: contigRecVector) {
                
                // std::cout << " CONTIG: " << contig.seq << std::endl;
                int match = 2;
                int mismatch = -6;
                int gap_open = -18;
                int gap_extend = 0;
                int match_special = 0;
                int mismatch_special = -20;

                AlignmentResult aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
                // AlignmentResult aln = align_with_ssw(contig.seq, regionalFasta);
                // std::cout << aln.query_start << " " << aln.query_end << " " << aln.ref_start << " " << aln.ref_end << std::endl;
                // std::cout << aln2.query_start << " " << aln2.query_end << " " << aln2.ref_start << " " << aln2.ref_end << std::endl;
                // std::cout << aln.cigarExtended << std::endl;
                // std::cout << aln2.cigarExtended << std::endl << std::endl;

                std::string alnCigar = aln.cigarExtended;
                int numOperations = aln.non_match_operations;
                int alnLength = aln.ref_end-aln.ref_start;
                double percLength = (double)alnLength/contig.seq.size();
                if (percLength < 0.9) {
                    continue;
                }

                bool hasOP = hasOpCluster(alnCigar, 3);
                if (hasOP) {
                    match_special = -10;
                    aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
                    // aln = aln2;
                    alnCigar = aln.cigarExtended;
                    numOperations = aln.non_match_operations;
                }

                std::vector<MutationSegment> mutationSegments;
                bool flag = false;
                int start = -1;
                int end = -1;
                const int maxDistance = 20;

                int numOp = 0;
                char prevOp;
                int mod = 0;

                bool isvalid = has_valid_end_operations(aln.cigarExtended);
                if (!isvalid) {
                    continue;
                }
                for (size_t i = 0; i < aln.cigarExtended.length(); ++i) {
                    char op = aln.cigarExtended[i];
                    if (!flag && (op == 'X' || op == 'D' || op == 'I')) {
                        flag = true;
                        start = i - mod;
                        end = i - mod;
                        numOp = 1;
                        prevOp = op;
                        mod = 0;
                    } 
                    else if (flag) {
                        if (op == 'X' || op == 'D' || op == 'I') {
                            if (i - end <= maxDistance) {
                                end = i;
                                mod = 0;
                            } else {
                                if (numOp > 1 || (numOp == 1 && prevOp != 'X')) {
                                    mutationSegments.push_back({start, end});
                                }
                                start = i - mod;
                                end = i - mod;
                                numOp = 1;
                                flag = true;
                                prevOp = op;
                                mod = 0;
                            }
                            numOp++;
                        } else {
                            if (numOp > 1 || (numOp == 1 && prevOp != 'X')) {
                                mutationSegments.push_back({start, end});
                            }
                            flag = false;
                        }
                    }
                }
                if (flag && (numOp > 1 || (numOp == 1 && prevOp != 'X'))) {
                    mutationSegments.push_back({start, end});
                }

                for (const auto& seg : mutationSegments) {
                    int newstart = seg.start;
                    int newend = seg.end;
                    std::string refAllele;
                    std::string altAllele;

                    bool isOutOfBound = false;
                    if (newstart < 0 || numOperations > 6 || ((newend - newstart) < 2 && numOperations > 6)) {
                        continue;
                    }

                    std::string reducedCigar = alnCigar.substr(newstart - 1, newend + 1 - newstart + 1);
                    // std::cout << reducedCigar << std::endl;
                    int i = aln.query_start + newstart;
                    // int j = newstart;
                    int j = newstart - 1;

                    for (auto& op : reducedCigar) {
                        char ref_ntd = regionalFasta[i];
                        char alt_ntd = aln.seq2_align[j];
                        if (i >=regionalFasta.size()-1) {
                            isOutOfBound = true;
                        }
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
                            } else if (op == 'I') {
                                if (alt_ntd != '-') {
                                    altAllele += alt_ntd;
                                }
                                j++;
                            }
                        }
                    }
                    if (refAllele == altAllele) {
                        continue;
                    }
                    if (isOutOfBound) {
                        continue;
                    }
                    // std::cout <<"alt:"<< altAllele << std::endl;
                    variant_t variant;
                    long newPos = aln.query_start + refStart + newstart - flankSpace;

                    variant.chr = contig.chr;
                    variant.pos = newPos;
                    variant.ref = refAllele;
                    variant.alt = altAllele;


                    // Initialize new readSupport count based on how many reads support the alt allele
                    int accurateReadSupport = 0;
                    std::unordered_set<std::string> usedReads2;
                    int plusStrand = 0;
                    int minusStrand = 0;

                    // contig.contigReads
                    std::mutex mtx;  // Declare a mutex for synchronizing access to shared resources

                    for (int i = 0; i < contig.contigReads.size(); i++) {
                        std::string readName = contig.contigReads[i].readName;  // Retrieve read name
                        std::string readSeq = contig.contigReads[i].read.Seq(); // Retrieve the corrected read sequence
                        std::string tag = readSeq + " " + readName;
                        if (threadLocalUsedReads2.find(tag) != threadLocalUsedReads2.end()) {
                            continue;  // Skip if already used in this thread
                        }

                        // {
                        //     // Lock the mutex before accessing the shared usedReads2 set
                        //     std::lock_guard<std::mutex> lock(mtx);

                        //     // Skip this read if it has already been used
                        //     if (usedReads2.find(tag) != usedReads2.end()) {
                        //         continue;
                        //     }
                        // }

                        // Get the start and end position of the read
                        long readStart = contig.contigReads[i].pos;              // Assuming Position() returns the read start position
                        long readEnd = readStart + readSeq.length();             // Read end is calculated based on sequence length

                        // Check if the read covers the variant position
                        if (!(readStart <= variant.pos && readEnd >= variant.pos)) {
                            continue;  // Skip if the read does not overlap the variant position
                        }

                        // Lock the mutex when updating shared variables
                        // {
                            std::lock_guard<std::mutex> lock(mtx);

                            // If the alternate allele is present in the read, count it toward read support
                            accurateReadSupport++;
                            // usedReads2.insert(tag);  // Mark this read as used
                            threadLocalUsedReads2.insert(tag);

                            // Update the strand-specific support counts
                            if (contig.contigReads[i].strand == '+') {
                                plusStrand++;
                            } else {
                                minusStrand++;
                            }
                        // }
                    }



                    variant.readSupport = accurateReadSupport;
                    variant.status = '.';
                    std::cout << variant.chr << " " <<  variant.pos << " " <<  variant.ref << " " << variant.alt << " " << variant.readSupport << std::endl;

                    AlignmentStats alnStats = CalculateAlignmentStats(bamFile, variant.chr, variant.pos);

                    variant.errorRate = alnStats.errorRate;
                    variant.chimericRate = alnStats.chimericRate;
                    variant.softClippedRate = alnStats.softClippedRate;

                    variant.mapQual = calculateMean(mapQualVector);

                    variant.totalDepth =alnStats.medianCoverage;
                    // variant.totalDepth = CalculateTotalDepth(bamFile, variant.chr, variant.pos);
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

                    // variant.plusStrand = plusStrand;
                    // variant.minusStrand = minusStrand;
                    // variant.strandBias = ComputeStrandBias(plusStrand, minusStrand);

                    uniqueVariants.push(variant);
                }
            }

            // if (contigsOut) {
            //     locusData.push_back(newLocus);
            // }
           
        }
    }

    std::vector<variant_t> result;
    variant_t var;
    while (uniqueVariants.try_pop(var)) {
        result.push_back(var);
    }

    return result;
}

