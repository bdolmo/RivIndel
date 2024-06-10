#include <iostream>
#include <string>
#include <fstream> 
#include <htslib/sam.h>
#include "BamRecord.h"
#include "BamReader.h"
#include "BamWriter.h"
#include <vector>
#include <map>
#include <optional>
#include "cluster_t.h"
#include <regex>
#include <filesystem>

struct TargetRegion {
    std::string chromosome;
    int start;
    int end;

    bool operator<(const TargetRegion& other) const {
        if (chromosome == other.chromosome) {
            return start < other.start;
        }
        return chromosome < other.chromosome;
    }
};

std::vector<TargetRegion> parseBedFile(const std::string& targetsBed) {
    std::vector<TargetRegion> targetRegions;
    std::ifstream bedFile(targetsBed);
    std::string line;
    while (std::getline(bedFile, line)) {
        std::istringstream iss(line);
        std::string chromosome;
        int start, end;
        if (!(iss >> chromosome >> start >> end)) {
            continue;
        }
        targetRegions.push_back({chromosome, start, end});
    }
    std::sort(targetRegions.begin(), targetRegions.end());
    return targetRegions;
}


std::vector<BamRecord> GetSoftClippedReadsInRegion(BamReader& reader, const std::string& region, int minSoftClip = 10) {
    std::vector<BamRecord> softClippedReads;
    reader.SetRegion(region);
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
                    int length = 0;
                    while (isdigit(c)) {
                        length = length * 10 + (c - '0');
                        c = cigar[++i];
                    }
                    if (c == 'S') {
                        if (leading) {
                            leadingSoftClips = length;
                        } else {
                            trailingSoftClips = length;
                        }
                    } else {
                        leading = false;
                    }
                }
            }

            if (leadingSoftClips >= minSoftClip || trailingSoftClips >= minSoftClip) {
                softClippedReads.push_back(record);
            }
        }
    }
    return softClippedReads;
}

// std::multimap<std::string, TargetRegion> parseBedFile(const std::string& targetsBed) {
//     std::multimap<std::string, TargetRegion> targetRegions;
//     std::ifstream bedFile(targetsBed);  // Ensure <fstream> is included
//     std::string line;
//     while (std::getline(bedFile, line)) {
//         std::istringstream iss(line);
//         std::string chromosome;
//         int start, end;
//         if (!(iss >> chromosome >> start >> end)) {
//             continue;
//         }
//         targetRegions.insert({chromosome, {start, end}});
//     }
//     return targetRegions;
// }

bool isOverlap(int start1, int end1, int start2, int end2) {
    return (start1 <= end2 && end1 >= start2);
}

bool isReadInTargetRegion(const BamRecord& record, const std::vector<TargetRegion>& targetRegions) {
    auto cmp = [](const TargetRegion& lhs, const TargetRegion& rhs) {
        return lhs.chromosome == rhs.chromosome ? lhs.start < rhs.start : lhs.chromosome < rhs.chromosome;
    };

    TargetRegion queryRegion = {record.chrName(), record.Position(), record.Position() + record.Seq().length()};
    auto lower = std::lower_bound(targetRegions.begin(), targetRegions.end(), queryRegion, cmp);

    // Check for overlap within a small window around the found region
    for (auto it = lower; it != targetRegions.end() && it->chromosome == record.chrName(); ++it) {
        if (isOverlap(record.Position(), record.Position() + record.Seq().length(), it->start, it->end)) {
            return true;
        }
        if (it->start > record.Position() + record.Seq().length()) {
            break; // No need to check further
        }
    }

    return false;
}

bool isValidReadSignal(BamRecord& record) {
    std::string cigar = record.cigarString();
    std::string mdTag = record.MDtag();

    // Counters and flags for CIGAR string analysis
    int lastOpPos = 0;
    int currentOpLength = 0;
    char lastOpType = ' ';
    bool foundIndel = false;
    bool foundTwoCloseIndels = false;

    int minOpSize = 1;
    int maxSeparation = 10;

    // Parsing CIGAR string
    for (char ch : cigar) {
        if (isdigit(ch)) {
            currentOpLength = currentOpLength * 10 + (ch - '0');
        } 
        else {
            if ((ch == 'I' || ch == 'D') && currentOpLength >= minOpSize) {
                if (foundIndel && (lastOpPos - currentOpLength <= maxSeparation)) {
                    foundTwoCloseIndels = true;
                } 
                else {
                    foundIndel = true;
                    lastOpPos = 0;
                }
            }
            lastOpType = ch;
            currentOpLength = 0;
        }
        lastOpPos += (lastOpType == 'M' || lastOpType == 'D' || lastOpType == 'I') ? currentOpLength : 0;
    }

    // Counter for MD tag analysis, considering mismatches only
    int mismatchCount = 0;
    for (size_t i = 0; i < mdTag.size(); ++i) {
        char ch = mdTag[i];
        if (!isdigit(ch) && ch != '^') {
            mismatchCount++;
            // If a mismatch is found within 5 bases after an indel
            if (foundIndel && lastOpPos <= 5) {
                return true;
            }
        }
        else if (isdigit(ch)) {
            int num = 0;
            while (i < mdTag.size() && isdigit(mdTag[i])) {
                num = num * 10 + (mdTag[i] - '0');
                ++i;
            }
            --i; // Compensate for the outer loop increment
            lastOpPos += num;
        }
    }

    // Skip reads with only a single mismatch
    if (mismatchCount == 1) {
        return false;
    }

    return foundIndel;
    // return foundTwoCloseIndels;
}


std::map<std::string, std::vector<clustered_aln_t>> clusterInformativeReads(
    const std::string bamFile, int minOpSize, int maxSeparation, int minSupport, 
    const std::string& tempBamFile, const std::string& targetsBed) {

    BamReader reader(bamFile);
    std::map<std::string, std::vector<clustered_aln_t>> candidateClusters;
    BamRecord record;
    int totalRecords = 0;
    std::string lastChromosome = "";

    auto targetRegions = parseBedFile(targetsBed);

    // Check if the temporary BAM file exists
    if (std::filesystem::exists(tempBamFile)) {
        std::cout << " INFO: Using existing temporary BAM file: " << tempBamFile << std::endl;
        BamReader tempReader(tempBamFile);
        while (tempReader.GetNextRecord(record)) {
            if (!isReadInTargetRegion(record, targetRegions)) {
                continue;
            }
            int windowIndex = record.Position() / 100;
            std::string key = record.chrName() + ":" + std::to_string(windowIndex);
            clustered_aln_t clusteredRead;
            clusteredRead.chromosome = record.chrName();
            clusteredRead.pos = record.Position();
            clusteredRead.read = record;
            clusteredRead.strand = record.GetStrand();

            double meanBaseQual = record.MeanBaseQuality();
            clusteredRead.mean_bq = meanBaseQual;

            candidateClusters[key].push_back(clusteredRead);
        }
    } else {
        std::cout << " INFO: Creating temporary BAM file: " << tempBamFile << std::endl;
        BamWriter tempWriter(tempBamFile, reader.getHeader());

        while (reader.GetNextRecord(record)) {
            if (!isReadInTargetRegion(record, targetRegions)) {
                continue;
            }
            totalRecords++;
            if (record.chrName() != lastChromosome && record.chrName() != "*") {
                std::cout << " INFO: Parsing " << record.chrName() << std::endl;
                lastChromosome = record.chrName();
            }

            if (isValidReadSignal(record)) {
                tempWriter.WriteRecord(record); // Write the valid read to the temporary BAM file

                int windowIndex = record.Position() / 100;
                std::string key = record.chrName() + ":" + std::to_string(windowIndex);
                clustered_aln_t clusteredRead;
                clusteredRead.chromosome = record.chrName();
                clusteredRead.pos = record.Position();
                clusteredRead.read = record;
                clusteredRead.strand = record.GetStrand();

                double meanBaseQual = record.MeanBaseQuality();
                clusteredRead.mean_bq = meanBaseQual;

                candidateClusters[key].push_back(clusteredRead);
            }
        }
        tempWriter.CreateIndex();
    }

    // Remove clusters that do not meet the minimum read support and add nearby soft-clipped reads
    BamReader softClipReader(bamFile);
    for (auto it = candidateClusters.begin(); it != candidateClusters.end();) {
        if (it->second.size() < minSupport) {
            it = candidateClusters.erase(it);
        } else {
            std::string chr = it->second[0].chromosome;
            int start = it->second[0].pos - 100;
            int end = it->second[0].pos + 100;
            std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);

            auto softClippedReads = GetSoftClippedReadsInRegion(softClipReader, region);

            for (const auto& scRead : softClippedReads) {
                clustered_aln_t clusteredRead;
                clusteredRead.chromosome = scRead.chrName();
                clusteredRead.pos = scRead.Position();
                clusteredRead.read = scRead;
                clusteredRead.strand = scRead.GetStrand();

                double meanBaseQual = scRead.MeanBaseQuality();
                clusteredRead.mean_bq = meanBaseQual;

                it->second.push_back(clusteredRead);
            }
            ++it;
        }
    }

    // return candidateClusters;
    // // Remove clusters that do not meet the minimum read support
    // for (auto it = candidateClusters.begin(); it != candidateClusters.end();) {
    //     if (it->second.size() < minSupport) {
    //         it = candidateClusters.erase(it);
    //     } else {
    //         ++it;
    //     }
    // }

    return candidateClusters;
}

std::string createTempBamFilename(const std::string& bamFile) {
    std::filesystem::path bamPath(bamFile);
    std::string baseName = bamPath.stem().string(); // Get the basename without extension
    return bamPath.parent_path().string() + "/" + baseName + ".tmp.bam"; // Append .tmp.bam suffix
}

// std::map<std::string, std::vector<clustered_aln_t>> clusterInformativeReads(
//     const std::string bamFile, int minOpSize, int maxSeparation, int minSupport, 
//     const std::string& tempBamFile, const std::string& targetsBed) {

//     BamReader reader(bamFile);
//     std::map<std::string, std::vector<clustered_aln_t>> candidateClusters;
//     BamRecord record;
//     int totalRecords = 0;
//     std::string lastChromosome = "";

//     std::vector<TargetRegion> targetRegions;
//     bool useTargets = !targetsBed.empty();
//     if (useTargets) {
//         targetRegions = parseBedFile(targetsBed);
//     }
    
//     // Check if the temporary BAM file exists
//     if (std::filesystem::exists(tempBamFile)) {
//         std::cout << " INFO: Using existing temporary BAM file: " << tempBamFile << std::endl;
//         BamReader tempReader(tempBamFile);
//         if (useTargets) {
//             for (const auto& region : targetRegions) {
//                 std::string regionStr = region.chromosome + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
//                 tempReader.SetRegion(regionStr);

//                 std::cout << "Region: " << regionStr << std::endl;

//                 while (tempReader.GetNextRecord(record)) {
//                     int windowIndex = record.Position() / 100;
//                     std::string key = record.chrName() + ":" + std::to_string(windowIndex);
//                     clustered_aln_t clusteredRead;
//                     clusteredRead.chromosome = record.chrName();
//                     clusteredRead.pos = record.Position();
//                     clusteredRead.read = record;
//                     clusteredRead.strand = record.GetStrand();
//                     double meanBaseQual = record.MeanBaseQuality();
//                     clusteredRead.mean_bq = meanBaseQual;
//                     candidateClusters[key].push_back(clusteredRead);
//                 }
//             }
//         }
//         else {

//             while (tempReader.GetNextRecord(record)) {
//                 int windowIndex = record.Position() / 100;
//                 std::string key = record.chrName() + ":" + std::to_string(windowIndex);
//                 clustered_aln_t clusteredRead;
//                 clusteredRead.chromosome = record.chrName();
//                 clusteredRead.pos = record.Position();
//                 clusteredRead.read = record;
//                 clusteredRead.strand = record.GetStrand();

//                 double meanBaseQual = record.MeanBaseQuality();
//                 clusteredRead.mean_bq = meanBaseQual;

//                 candidateClusters[key].push_back(clusteredRead);
//             }
//         }
//     } else {
//         std::cout << " INFO: Creating temporary BAM file: " << tempBamFile << std::endl;
//         BamWriter tempWriter(tempBamFile, reader.getHeader());

//         if (useTargets) {
//             for (const auto& region : targetRegions) {
//                 std::string regionStr = region.chromosome + ":" + std::to_string(region.start) + "-" + std::to_string(region.end);
//                 reader.SetRegion(regionStr);
//                 std::cout << "Region: " << regionStr << std::endl;

//                 while (reader.GetNextRecord(record)) {
//                     totalRecords++;
//                     if (record.chrName() != lastChromosome && record.chrName() != "*") {
//                         std::cout << " INFO: Parsing " << record.chrName() << std::endl;
//                         lastChromosome = record.chrName();
//                     }
//                     if (isValidReadSignal(record)) {
//                         tempWriter.WriteRecord(record); // Write the valid read to the temporary BAM file

//                         int windowIndex = record.Position() / 100;
//                         std::string key = record.chrName() + ":" + std::to_string(windowIndex);
//                         clustered_aln_t clusteredRead;
//                         clusteredRead.chromosome = record.chrName();
//                         clusteredRead.pos = record.Position();
//                         clusteredRead.read = record;
//                         clusteredRead.strand = record.GetStrand();
//                         double meanBaseQual = record.MeanBaseQuality();
//                         clusteredRead.mean_bq = meanBaseQual;
//                         candidateClusters[key].push_back(clusteredRead);
//                     }
//                 }
//             }
//         } 
//         else {
//             while (reader.GetNextRecord(record)) {
//                 totalRecords++;
//                 if (record.chrName() != lastChromosome && record.chrName() != "*") {
//                     std::cout << " INFO: Parsing " << record.chrName() << std::endl;
//                     lastChromosome = record.chrName();
//                 }

//                 if (isValidReadSignal(record)) {
//                     tempWriter.WriteRecord(record); // Write the valid read to the temporary BAM file

//                     int windowIndex = record.Position() / 100;
//                     std::string key = record.chrName() + ":" + std::to_string(windowIndex);
//                     clustered_aln_t clusteredRead;
//                     clusteredRead.chromosome = record.chrName();
//                     clusteredRead.pos = record.Position();
//                     clusteredRead.read = record;
//                     clusteredRead.strand = record.GetStrand();

//                     double meanBaseQual = record.MeanBaseQuality();
//                     clusteredRead.mean_bq = meanBaseQual;

//                     candidateClusters[key].push_back(clusteredRead);
//                 }
//             }
//         }

//         // tempWriter.CreateIndex();
//     }

//     // Remove clusters that do not meet the minimum read support
//     for (auto it = candidateClusters.begin(); it != candidateClusters.end();) {
//         if (it->second.size() < minSupport) {
//             it = candidateClusters.erase(it);
//         } else {
//             ++it;
//         }
//     }

//     return candidateClusters;
// }
