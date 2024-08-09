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
#include <set>

struct TargetRegion {
    std::string chromosome;
    int start;
    int end;

    // bool operator<(const TargetRegion& other) const {
    //     if (chromosome == other.chromosome) {
    //         return start < other.start;
    //     }
    //     return chromosome < other.chromosome;
    // }
	bool operator<(const TargetRegion& x) const {
		return (start < x.start);
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


std::map<std::string, std::vector<TargetRegion>> parseBedFile2(const std::string& targetsBed) {
    std::map<std::string, std::vector<TargetRegion>> targetRegionsMap;
    std::ifstream bedFile(targetsBed);
    std::string line;

    while (std::getline(bedFile, line)) {
        std::istringstream iss(line);
        std::string chromosome;
        int start, end;
        if (!(iss >> chromosome >> start >> end)) {
            continue;
        }
        targetRegionsMap[chromosome].push_back({chromosome, start, end});
    }

    // Sort the target regions for each chromosome
    for (auto& pair : targetRegionsMap) {
        std::sort(pair.second.begin(), pair.second.end());
    }

    return targetRegionsMap;
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
                // std::cout << cigar << " " << leadingSoftClips << " " << trailingSoftClips << std::endl;
                softClippedReads.push_back(record);
            }
        }
    }
    return softClippedReads;
}

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

    if (record.mapQual() < 20) {
        return false;
    }

    for (char ch : cigar) {
        if (isdigit(ch)) {
            currentOpLength = currentOpLength * 10 + (ch - '0');
        } 
        else {
            if ((ch == 'I' || ch == 'D') && currentOpLength >= minOpSize) {
                // Check if we found a single indel
                if (currentOpLength >= 1) {
                    // std::cout << cigar << " " << ch << " " << currentOpLength << std::endl;
                    return true;  // Found a single I or D operation
                }
            }
            lastOpType = ch;
            currentOpLength = 0;

        }
    }

    currentOpLength = 0;
    lastOpPos = 0;
    lastOpType = ' ';

    // Parsing CIGAR string
    for (char ch : cigar) {
        if (isdigit(ch)) {
            currentOpLength = currentOpLength * 10 + (ch - '0');
        }
        else {
            if ((ch == 'I' || ch == 'D') && currentOpLength >= minOpSize) {
                // Check if we found a single indel
                if (currentOpLength >= 1) {
                    return true;  // Found a single I or D operation
                }

                // Check if it's a valid indel and not too close to the previous one
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

bool containsN(const std::string& read) {
    return read.find('N') != std::string::npos;
}

bool binarySearch(const std::vector<TargetRegion>& vec, int n, int64_t start, int64_t end) {
    int low = 0;
    int high = n - 1;

    while (low <= high) {
        int mid = (low + high) / 2;
        const TargetRegion& region = vec[mid];

        // Check if the query range is within the current target region
        if (start >= region.start && end <= region.end) {
            return true;
        }

        // Narrow the search range
        if (end < region.start) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    return false;
}


bool liesOnTargetRegion(std::vector<TargetRegion>& badRegions, std::string& readChrom, int64_t readPos, int& readLen) {
    int64_t readEnd = readPos + readLen;
	bool doesOverlap = binarySearch( badRegions, badRegions.size(), readPos, readEnd );
	if (!doesOverlap) {
		return 0;
	}
	else {
		return 1;
	}
}


bool liesOnBlackRegion(std::vector<TargetRegion>& badRegions, std::string& readChrom, int64_t readPos, int& readLen) {
    int64_t readEnd = readPos + readLen;
	bool doesOverlap = binarySearch( badRegions, badRegions.size(), readPos, readEnd );
	if (!doesOverlap) {
		return 0;
	}
	else {
		return 1;
	}
}


void createTempBamFile(const std::string& bamFile, const std::string& chromosomeName, 
                       const std::string& tempBamFile, const std::string& targetsBed,
                       const std::string& excludeBed, int minOpSize, int maxSeparation, int minSupport) {

    BamReader reader(bamFile);
    reader.SetRegion(chromosomeName);  // Set the region to the specific chromosome
    BamRecord record;
    int totalRecords = 0;

    auto excludeRegions = parseBedFile2(excludeBed);
    auto targetRegions = parseBedFile2(targetsBed);

    std::cout << " INFO: Creating temporary BAM file: " << tempBamFile << std::endl;
    BamWriter tempWriter(tempBamFile, reader.getHeader());

    while (reader.GetNextRecord(record)) {
        int readLength = record.Seq().length();
        std::string chromosomeName = record.chrName();
        
        if (containsN(record.Seq())) {
            continue;
        }

        totalRecords++;

        if (isValidReadSignal(record)) {
            tempWriter.WriteRawRecord(record); // Write the valid read to the temporary BAM file
        }
    }
    tempWriter.Close();

    tempWriter.CreateIndex();

    // tempWriter.CreateIndex();
}

std::map<std::string, std::vector<clustered_aln_t>> clusterInformativeReads(
    const std::string& bamFile,
    const std::string& tempBamFile, int minOpSize, int maxSeparation, int minSupport, 
    const std::string& targetsBed, const std::string& excludeBed, const std::string& clustersLog, const std::string& chromosomeName) {

    BamReader reader(tempBamFile);
    std::map<std::string, std::vector<clustered_aln_t>> candidateClusters;
    BamRecord record;
    int totalRecords = 0;

    auto targetRegions = parseBedFile2(targetsBed);
    auto excludeRegions = parseBedFile2(excludeBed);

    std::ofstream clustersFile(clustersLog);

    bool inCluster = false;
    std::vector<clustered_aln_t> currentCluster;
    int clusterStart = 0;
    int clusterEnd = 0;

    std::set<std::string> processedReads;

while (reader.GetNextRecord(record)) {
    int readLength = record.Seq().length();
    std::string readChrName = record.chrName();
    std::string chromosomeName = record.chrName();
    std::string readName = record.Qname();

    // if (processedReads.find(readName) != processedReads.end()) {
    //     continue; // Skip if the read has already been processed
    // }

    if (readChrName != chromosomeName) {
        continue;
    }
    if (liesOnBlackRegion(excludeRegions[readChrName], chromosomeName, record.Position(), readLength)) {
        continue;
    }
    if (!liesOnTargetRegion(targetRegions[readChrName], chromosomeName, record.Position(), readLength)) {
        continue;
    }
    if (containsN(record.Seq())) {
        continue;
    }

    // processedReads.insert(readName); // Mark the read as processed

    int readPos = record.Position();
    int windowIndex = readPos / 100;
    std::string key = readChrName + ":" + std::to_string(windowIndex);

    clustered_aln_t clusteredRead;
    clusteredRead.chromosome = readChrName;
    clusteredRead.pos = readPos;
    clusteredRead.read = record;
    clusteredRead.strand = record.GetStrand();
    clusteredRead.mean_bq = record.MeanBaseQuality();

    if (!inCluster) {
        // Start a new cluster
        clusterStart = readPos - 100; // 500 bp upstream
        clusterEnd = readPos + 100;   // 500 bp downstream
        currentCluster.push_back(clusteredRead);
        inCluster = true;
    } else {
        // Add read to the current cluster if it's within the cluster range and the total cluster size is within 500 bp
        if (readPos >= clusterStart && readPos <= clusterEnd && (clusterEnd - clusterStart <= 500)) {
            currentCluster.push_back(clusteredRead);
            // Update the cluster end if this read extends beyond it
            clusterEnd = std::max(clusterEnd, readPos + 100);
        } else {
            // End the current cluster and start a new one
            if (currentCluster.size() >= minSupport) {
                candidateClusters[key] = currentCluster;
            }
            currentCluster.clear();
            clusterStart = readPos - 100;
            clusterEnd = readPos + 100;
            currentCluster.push_back(clusteredRead);
        }
    }
}

    // while (reader.GetNextRecord(record)) {
    //     int readLength = record.Seq().length();
    //     std::string readChrName = record.chrName();
    //     std::string chromosomeName = record.chrName();

    //     if (readChrName != chromosomeName) {
    //         continue;
    //     }
    //     if (liesOnBlackRegion(excludeRegions[readChrName], chromosomeName, record.Position(), readLength)) {
    //         continue;
    //     }
    //     if (!liesOnTargetRegion(targetRegions[readChrName], chromosomeName, record.Position(), readLength)) {
    //         continue;
    //     }
    //     if (containsN(record.Seq())) {
    //         continue;
    //     }
    //     // std::string readName = record.Qname() + record.Seq();

    //     // if (processedReads.find(readName) != processedReads.end()) {
    //     //     continue; // Skip if the read has already been processed
    //     // }

    //     // processedReads.insert(readName); // Mark the read as processed

    //     int readPos = record.Position();
    //     int windowIndex = readPos / 100;
    //     std::string key = readChrName + ":" + std::to_string(windowIndex);
        
    //     clustered_aln_t clusteredRead;
    //     clusteredRead.chromosome = readChrName;
    //     clusteredRead.pos = readPos;
    //     clusteredRead.read = record;
    //     clusteredRead.strand = record.GetStrand();
    //     clusteredRead.mean_bq = record.MeanBaseQuality();

    //     if (!inCluster) {
    //         // Start a new cluster
    //         clusterStart = readPos - 100; // 500 bp upstream
    //         clusterEnd = readPos + 100;   // 500 bp downstream
    //         currentCluster.push_back(clusteredRead);
    //         inCluster = true;
    //     } else {
    //         // Add read to the current cluster if it's within the cluster range
    //         if (readPos >= clusterStart && readPos <= clusterEnd) {
    //             currentCluster.push_back(clusteredRead);
    //             // Update the cluster end if this read extends beyond it
    //             clusterEnd = std::max(clusterEnd, readPos + 100);
    //         } else {
    //             // End the current cluster and start a new one
    //             if (currentCluster.size() >= minSupport) {
    //                 candidateClusters[key] = currentCluster;
    //             }
    //             currentCluster.clear();
    //             clusterStart = readPos - 100;
    //             clusterEnd = readPos + 100;
    //             currentCluster.push_back(clusteredRead);
    //         }
    //     }
        
    // }

    // Add soft-clipped reads to the clusters
    std::cout << " INFO: Rescueing soft-clipped reads for " << chromosomeName << std::endl;

    // Remove clusters that do not meet the minimum read support and add nearby soft-clipped reads
    BamReader softClipReader(bamFile);
    for (auto it = candidateClusters.begin(); it != candidateClusters.end();) {
        if (it->second.size() < minSupport) {
            it = candidateClusters.erase(it);
        }
        else {
            int64_t minPos = 9000000000000000;
            int64_t maxPos = 0;
            for (auto& i : it->second) {
                int64_t rStart = i.pos;
                int64_t rEnd = i.pos+i.read.Seq().length();

                if (rStart <=  minPos) {
                    minPos = rStart;
                }
                if (rEnd >= maxPos+i.read.Seq().length()) {
                    maxPos = rEnd;
                }
            }

            std::string chr = it->second[0].chromosome;
            std::string region = chr + ":" + std::to_string(minPos) + "-" + std::to_string(maxPos);

            auto softClippedReads = GetSoftClippedReadsInRegion(softClipReader, region);
            for (const auto& scRead : softClippedReads) {
                clustered_aln_t clusteredRead;
                clusteredRead.chromosome = scRead.chrName();
                clusteredRead.pos = scRead.Position();
                clusteredRead.read = scRead;
                clusteredRead.strand = scRead.GetStrand();
                std::string readName = scRead.Qname() + scRead.Seq();

                // if (processedReads.find(readName) != processedReads.end()) {
                //     continue; // Skip if the read has already been processed
                // }

                double meanBaseQual = scRead.MeanBaseQuality();
                clusteredRead.mean_bq = meanBaseQual;
                // if (processedReads.find(readName) == processedReads.end()) {
                //     it->second.push_back(clusteredRead);
                //     processedReads.insert(readName);
                // }
                it->second.push_back(clusteredRead);
            }
 
            clustersFile << chr << "\t" << minPos << "\t" << maxPos << "\t" << it->second.size() << "\n";
            ++it;
        }
    }
    clustersFile.close();

    return candidateClusters;
}