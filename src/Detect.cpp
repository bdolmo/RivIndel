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
        std::vector<std::string> kmers = extractKmers(sequence, k);
        for (const auto& kmer : kmers) {
            kmerCounts[kmer]++;
        }
    }
    return kmerCounts;
}

// Function to identify and correct erroneous k-mers
std::string correctSequence(const std::string& sequence, int k, const std::unordered_map<std::string, int>& kmerCounts, int threshold) {
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
        correctedSequences.push_back(correctSequence(sequence, k, kmerCounts, threshold));
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
    EdlibAlignResult result = edlibAlign(read.c_str(), read.size(), contig.c_str(), contig.size(),
        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

    int editDistance = result.status == EDLIB_STATUS_OK ? result.editDistance : -1;
    std::string cigar;

    if (result.status == EDLIB_STATUS_OK && result.alignment) {
        char* cigarC = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        cigar = std::string(cigarC);
        free(cigarC);
    }


    // Print out the offsets and the CIGAR string
    if (result.status == EDLIB_STATUS_OK) {
        // std::cout << "Read Start: " << result.startLocations[0] << std::endl;
        // std::cout << "Read End: " << result.endLocations[0] << std::endl;
        // std::cout << "Contig Start: " << result.startLocations[1] << std::endl;
        // std::cout << "Contig End: " << result.endLocations[1] << std::endl;
        // std::cout << read << " " << contig << std::endl;
    // printf("%d\n", result.editDistance);
    // printf("%d\n", result.alignmentLength);
    // printf("%d\n", result.endLocations[0]);
        // std::cout << result.alignmentLength << std::endl;

        // std::cout << "Editdistance: " << editDistance << " CIGAR String: " << cigar << " " << result.alignment << std::endl << std::endl;
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

struct AlignmentStats {
    double errorRate;
    double chimericRate;
    double softClippedRate;
};

// Function to calculate various alignment statistics
AlignmentStats CalculateAlignmentStats(const std::string& bamFile, const std::string& chr, int position) {
    int totalBases = 0;
    int mismatches = 0;
    int totalReads = 0;
    int chimericReads = 0;
    int softClippedReads = 0;
    BamReader newreader(bamFile); 

    std::string region = chr + ":" + std::to_string(position - 100) + "-" + std::to_string(position + 100);
    newreader.SetRegion(region);
    int maxReads = 500;
    int count = 0;
    BamRecord record;
    while (newreader.GetNextRecord(record)) {
        totalReads++;
        std::string mdTag = record.MDtag();
        totalBases += record.Seq().length();
        count++;

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
                --i; // Correct the position after exiting the while loop
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

    AlignmentStats stats;
    if (totalBases == 0) {
        stats.errorRate = 0.0;
    } else {
        stats.errorRate = static_cast<double>(mismatches) / totalBases;
    }
    stats.chimericRate = static_cast<double>(chimericReads) / totalReads;

    stats.softClippedRate = static_cast<double>(softClippedReads) / totalReads;

    return stats;
}


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

    // Optionally, set custom parameters for error correction
    // opt.ec_k = 11;        // Set k-mer length for error correction
    // opt.min_asm_ovlp = 11; // Set minimum overlap length during assembly
    // opt.min_merge_len = 11; // Set minimum length for merging overlaps

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
            std::cout << " INFO: Interrogating "<< chr << ":" << refStart << "-" << refEnd << std::endl;

            Locus newLocus;
            newLocus.chromosome = chr;
            newLocus.start = refStart;
            newLocus.end = refEnd;

            std::string regionalFasta = ref.fetchSequence(chr, refStart - 49, refEnd + 50);
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
                // std::cout << "total_reads: " << reads.size() << std::endl;
                for (int i=0; i<reads.size();i++) {
                    // std::string readSeq = reads[i].read.Seq();
                    std::string readSeq = corrReads[i];
                    // if (usedReads.find(readSeq) != usedReads.end()) {
                    //     // Skip this read if it has already been used
                    //     continue;
                    // }
                    mapQualVector.push_back(reads[i].read.mapQual());
                    auto [editDistance, cigar] = alignReadToContig(readSeq, contig);
                    if (editDistance >= 0 && editDistance <= 10) {
                        readSupport++;
                        isGood = true;
                        contigStart = std::min(contigStart, static_cast<long>(reads[i].pos));
                        contigEnd = std::max(contigEnd, static_cast<long>(reads[i].pos + readSeq.size()));
                        contigRecord.seq = contig;

                        if (reads[i].strand == '+') {
                           plusStrand++;
                        } else {
                           minusStrand++;
                        }
                        usedReads.insert(readSeq);
                    }
                    else{
                        std::string revContig = reverseComplement(contig);
                        auto [editDistance, cigar] = alignReadToContig(readSeq, revContig);
                        if (editDistance >= 0 && editDistance <= 10) {
                            readSupport++;
                            isGood = true;
                            contigStart = std::min(contigStart, static_cast<long>(reads[i].pos));
                            contigEnd = std::max(contigEnd, static_cast<long>(reads[i].pos + readSeq.size()));
                            contigRecord.seq = revContig;

                            if (reads[i].strand == '+') {
                            plusStrand++;
                            } else {
                            minusStrand++;
                            }
                            usedReads.insert(readSeq);
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
                int alnLength = aln.ref_end-aln.ref_start;
                double percLength = (double)alnLength/contig.seq.size();
                //  std::cout << contig.seq << " " << alnCigar  << " " << numOperations << " "  << percLength << std::endl;
                if (percLength < 0.9) {
                    continue;
                }

                bool hasOP = hasOpCluster(alnCigar, 3);
                if (hasOP) {
                    match_special = -10;
                    aln = affine_semiglobal(regionalFasta, contig.seq, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
                    alnCigar = aln.cigarExtended;
                    numOperations = aln.non_match_operations;
                }

                std::vector<MutationSegment> mutationSegments;
                // std::vector<MutationSegment> localMutationSegments;
                bool flag = false;
                int start = -1;
                int end = -1;
                const int maxDistance = 20;

                int numOp = 0;
                char prevOp;
                int mod = 0;

                for (size_t i = 0; i < aln.cigarExtended.length(); ++i) {
                    char op = aln.cigarExtended[i];
                    // std::cout << "operation " << op << std::endl;

                    if (!flag && (op == 'X' || op == 'D' || op == 'I')) {
                        flag = true;
                        start = i - mod;
                        end = i - mod;
                        numOp = 1;
                        prevOp = op;
                        mod = 0;
                    } else if (flag) {
                        // std::cout << "flag in " << i << " " << end  << " " << std::endl;

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

                    if (newstart < 0 || numOperations > 6 || ((newend - newstart) < 2 && numOperations > 6)) {
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

                    variant_t variant;
                    long newPos = aln.query_start + refStart + newstart - 49;

                    variant.chr = contig.chr;
                    variant.pos = newPos;
                    variant.ref = refAllele;
                    variant.alt = altAllele;
                    variant.readSupport = contig.readSupport;
                    variant.status = '.';

                    AlignmentStats alnStats = CalculateAlignmentStats(bamFile, variant.chr, variant.pos);

                    variant.errorRate = alnStats.errorRate;
                    variant.chimericRate = alnStats.chimericRate;
                    variant.softClippedRate = alnStats.softClippedRate;

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

