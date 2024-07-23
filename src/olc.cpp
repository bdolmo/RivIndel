#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Structure to represent an overlap
struct Overlap {
    string read;  // The overlapping read
    int length;   // The length of the overlap
};

// Function to calculate k-mer diversity
float calculateKmerDiversity2(const string& sequence, int k = 2) {
    unordered_set<string> uniqueKmers;
    int n = sequence.size();

    // Extract all k-mers and add to the set
    for (int i = 0; i <= n - k; ++i) {
        string kmer = sequence.substr(i, k);
        uniqueKmers.insert(kmer);
    }

    // Calculate diversity as the number of unique k-mers divided by the total possible k-mers
    int totalKmers = n - k + 1;
    if (totalKmers <= 0) return 0.0f; // Handle edge case where sequence length < k
    float diversity = static_cast<float>(uniqueKmers.size()) / totalKmers;

    return diversity;
}

// Function to find the best overlap between two sequences with a minimum overlap length, allowing up to K mismatches, and filtering low complexity sequences
int findBestOverlap(const string& seq1, const string& seq2, int minOverlap) {
    float minDiversity = 0.25;
    int maxMismatches = 1;
    int maxOverlap = 0;
    int len1 = seq1.length();
    int len2 = seq2.length();

    for (int i = minOverlap; i <= len1 && i <= len2; ++i) {
        string suffix = seq1.substr(len1 - i);
        string prefix = seq2.substr(0, i);

        // Check for exact match in the initial overlap
        if (suffix == prefix) {
            float suffixDiv = calculateKmerDiversity2(suffix);

            // Check for low complexity sequences
            if (suffixDiv < minDiversity) {
                break;
            }

            maxOverlap = i;
            int mismatches = 0;
            int j = i;

            // Extend comparison beyond the initial overlap region allowing up to maxMismatches
            while (j < len1 && j < len2) {
                if (seq1[len1 - i + j] != seq2[j]) {
                    mismatches++;
                    if (mismatches > maxMismatches) {
                        break;
                    }
                }
                j++;
                maxOverlap = j; // Update maxOverlap to the extended region
            }
        }
    }

    return maxOverlap;
}

// Function to build layout using a greedy approach
unordered_map<string, vector<Overlap>> buildGreedyLayout(const vector<string>& sequences, int minOverlap, int minAllowedOverlap) {
    int n = sequences.size();
    unordered_map<string, vector<Overlap>> overlapMap;  // Map to store overlaps
    vector<bool> used(n, false);     // Flag for used reads

    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;       // Skip used reads
        int current = i;
        used[current] = true;

        while (true) {
            int bestOverlap = 0;
            int bestNext = -1;

            for (int j = 0; j < n; ++j) {
                if (used[j] || current == j) continue;
                int overlap = findBestOverlap(sequences[current], sequences[j], minOverlap);
                if (overlap > bestOverlap && overlap >= minAllowedOverlap) {
                    bestOverlap = overlap;
                    bestNext = j;
                }
            }

            if (bestNext == -1) break; // No more overlaps found
            used[bestNext] = true;

            // Store the overlap in the map
            overlapMap[sequences[current]].push_back({sequences[bestNext], bestOverlap});
            current = bestNext;
        }
    }

    // Print the content of overlapMap
    // cout << "Overlap Map:" << endl;
    // for (const auto& from : overlapMap) {
    //     cout << "Read " << from.first << " has overlaps with:" << endl;
    //     for (const auto& to : from.second) {
    //         cout << "    Read " << to.read << " with overlap of " << to.length << " bases." << endl;
    //     }
    // }

    return overlapMap;
}

// Function to remove low-weight arcs
void removeLowWeightArcs(unordered_map<string, vector<Overlap>>& overlapMap, int weightThreshold) {
    for (auto& [read, overlaps] : overlapMap) {
        overlaps.erase(remove_if(overlaps.begin(), overlaps.end(), [weightThreshold](const Overlap& o) {
            return o.length < weightThreshold;
        }), overlaps.end());
    }
}

// Function to remove transitively inferable edges
void removeTransitivelyInferableEdges(unordered_map<string, vector<Overlap>>& overlapMap) {
    for (auto& [read, overlaps] : overlapMap) {
        unordered_set<string> directConnections;
        for (const auto& o : overlaps) {
            directConnections.insert(o.read);
        }

        vector<Overlap> toRemove;
        for (const auto& o1 : overlaps) {
            for (const auto& o2 : overlapMap[o1.read]) {
                if (directConnections.count(o2.read)) {
                    toRemove.push_back(o2);
                }
            }
        }

        for (const auto& tr : toRemove) {
            overlaps.erase(remove_if(overlaps.begin(), overlaps.end(), [&tr](const Overlap& o) {
                return o.read == tr.read;
            }), overlaps.end());
        }
    }
}

// Function to generate the consensus sequence from the layout for a single component
string generateConsensus(const vector<string>& sequences, const vector<int>& layout, unordered_map<string, vector<Overlap>>& overlapMap) {
    if (layout.empty()) return "";

    string consensus = sequences[layout[0]];
    for (size_t i = 1; i < layout.size(); ++i) {
        string prevRead = sequences[layout[i - 1]];
        string currRead = sequences[layout[i]];



        int overlap = 0;

        // Find the overlap length from the map
        for (const auto& o : overlapMap[prevRead]) {
            if (o.read == currRead) {
                overlap = o.length;
                break;
            }
        }

        std::cout << "PrevRead: " << prevRead << std::endl;
        std::cout << "currRead: " << currRead << std::endl;
        std::cout << "Overlap: " << overlap << std::endl;
        std::cout << "prevconsensus" << consensus<< std::endl;

        // Combine the overlapping regions
        for (int j = 0; j < overlap; ++j) {
            char base1 = consensus[consensus.size() - overlap + j];
            char base2 = currRead[j];
            if (base1 != base2) {
                // Take the most common nucleotide (simple voting)
                unordered_map<char, int> count;
                count[base1]++;
                count[base2]++;
                consensus[consensus.size() - overlap + j] = count[base1] >= count[base2] ? base1 : base2;
            }
        }

        consensus += currRead.substr(overlap-1);
              std::cout << consensus<< std::endl;

       std::cout << std::endl;


        // Append the non-overlapping part of the sequence
        
    }
    return consensus;
}

// Function to remove duplicate reads
vector<string> removeDuplicates(const vector<string>& reads) {
    unordered_set<string> uniqueReads(reads.begin(), reads.end());
    return vector<string>(uniqueReads.begin(), uniqueReads.end());
}

// Function to perform the OLC assembly
pair<vector<string>, vector<string>> OLC(const vector<string>& reads, int minOverlap, int weightThreshold, unordered_map<string, vector<Overlap>>& overlapMap) {
    int minAllowedOverlap = minOverlap;
    vector<string> uniqueReads = removeDuplicates(reads);
    vector<string> contigs;
    vector<string> unassembledReads;
    overlapMap = buildGreedyLayout(uniqueReads, minOverlap, minAllowedOverlap);

    // Remove low-weight arcs
    // removeLowWeightArcs(overlapMap, weightThreshold);

    // // Remove transitively inferable edges
    // removeTransitivelyInferableEdges(overlapMap);

    unordered_set<int> assembledReads;

    for (const auto& component : overlapMap) {
        vector<int> layout;
        auto it = find(uniqueReads.begin(), uniqueReads.end(), component.first);
        if (it != uniqueReads.end()) {
            int index = distance(uniqueReads.begin(), it);
            layout.push_back(index);
            assembledReads.insert(index);
        }
        for (const auto& overlap : component.second) {
            auto it = find(uniqueReads.begin(), uniqueReads.end(), overlap.read);
            if (it != uniqueReads.end()) {
                int index = distance(uniqueReads.begin(), it);
                layout.push_back(index);
                assembledReads.insert(index);
            }
        }
        string consensus = generateConsensus(uniqueReads, layout, overlapMap);
        contigs.push_back(consensus);
    }

    // Find unassembled reads
    for (int i = 0; i < uniqueReads.size(); ++i) {
        if (assembledReads.find(i) == assembledReads.end()) {
            unassembledReads.push_back(uniqueReads[i]);
        }
    }

    return {contigs, unassembledReads};
}

// int main() {
//     vector<string> sequences = {"AGCTAGCTAG", "CTAGCTAGCTA", "GCTAGCTAGCA", "TAGCTAGCTAG", "CTAGCTAGCTA"};
//     int minOverlap = 11;
//     int weightThreshold = 5; // Example threshold for removing low-weight arcs

//     unordered_map<string, vector<Overlap>> overlapMap;

//     // Assemble contigs
//     auto result = OLC(sequences, minOverlap, weightThreshold, overlapMap);
//     vector<string> contigs = result.first;
//     vector<string> unassembledReads = result.second;

//     // Print the contigs
//     cout << "Contigs:" << endl;
//     for (const auto& contig : contigs) {
//         cout << contig << endl;
//     }

//     // Print unassembled reads
//     cout << "Unassembled reads:" << endl;
//     for (const auto& read : unassembledReads) {
//         cout << read << endl;
//     }

//     // Print the overlap map
//     cout << "Overlap Map:" << endl;
//     for (const auto& from : overlapMap) {
//         for (const auto& to : from.second) {
//             cout << "Read " << from.first << " overlaps with Read " << to.read << " by " << to.length << " bases." << endl;
//         }
//     }

//     return 0;
// }
