#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>

using namespace std;

// Function to find the best overlap between two sequences with a minimum overlap length
int findBestOverlap(const string& seq1, const string& seq2, int minOverlap) {
    int maxOverlap = 0;
    int len1 = seq1.length();
    int len2 = seq2.length();

    for (int i = minOverlap; i <= len1 && i <= len2; ++i) {
        string suffix = seq1.substr(len1 - i);
        string prefix = seq2.substr(0, i);

        // Check for exact match
        if (suffix == prefix) {
            maxOverlap = i;
            // Extend comparison beyond the initial overlap region
            int j = i;
            while (j < len1 && j < len2 && seq1[len1 - i + j] == seq2[j]) {
                j++;
            }
            maxOverlap = j;
        }
    }

    return maxOverlap;
}

// Function to build layout using a greedy approach
vector<vector<int>> buildGreedyLayout(const vector<string>& sequences, int minOverlap, int minAllowedOverlap) {
    int n = sequences.size();
    vector<vector<int>> components;
    vector<bool> used(n, false);

    for (int i = 0; i < n; ++i) {
        if (used[i]) continue;
        used[i] = true;
        vector<int> component = {i};
        int current = i;
        bool foundOverlap = false;

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

            if (bestNext == -1) break;
            used[bestNext] = true;
            component.push_back(bestNext);
            foundOverlap = true;
            current = bestNext;
        }

        // Only add the component if it has valid overlaps
        if (foundOverlap) {
            components.push_back(component);
        }
    }

    return components;
}


#include <unordered_map>

// Function to generate the consensus sequence from the layout for a single component
string generateConsensus(const vector<string>& sequences, const vector<int>& layout, int minOverlap) {
    if (layout.empty()) return "";

    string consensus = sequences[layout[0]];
    for (size_t i = 1; i < layout.size(); ++i) {
        int overlap = findBestOverlap(sequences[layout[i - 1]], sequences[layout[i]], minOverlap);

        // Combine the overlapping regions
        for (int j = 0; j < overlap; ++j) {
            char base1 = consensus[consensus.size() - overlap + j];
            char base2 = sequences[layout[i]][j];
            if (base1 != base2) {
                // Take the most common nucleotide (simple voting)
                unordered_map<char, int> count;
                count[base1]++;
                count[base2]++;
                consensus[consensus.size() - overlap + j] = count[base1] >= count[base2] ? base1 : base2;
            }
        }

        // Append the non-overlapping part of the sequence
        consensus += sequences[layout[i]].substr(overlap);
    }
    return consensus;
}



pair<vector<string>, vector<string>> OLC(const vector<string>& reads, int minOverlap) {

    int minAllowedOverlap = minOverlap;
    vector<string> contigs;
    vector<string> unassembledReads;
    vector<vector<int>> components = buildGreedyLayout(reads, minOverlap, minAllowedOverlap);

    unordered_set<int> assembledReads;

    for (const auto& component : components) {
        if (component.size() > 1) {  // Ensure the component has more than one read
            string consensus = generateConsensus(reads, component, minOverlap);
            contigs.push_back(consensus);
            assembledReads.insert(component.begin(), component.end());
        } else {
            unassembledReads.push_back(reads[component[0]]);
        }
    }

    // Find unassembled reads
    for (int i = 0; i < reads.size(); ++i) {
        if (assembledReads.find(i) == assembledReads.end()) {
            unassembledReads.push_back(reads[i]);
        }
    }

    return {contigs, unassembledReads};
}

int main() {
    vector<string> sequences = {"ACGGGATATTACGTCGTCGATGCTAGCTAGTC", "GGATATTACGTCGTCGATGCTAGCTAGTCCCCCCC", "CGATGCTAGCTAGTCCCCCCCATTAACGATCGTA", "CCATTGGGGGGGGGGGGGGGGGGGGGG"};
    int minOverlap = 5;
    int minAllowedOverlap = 5;

    // Assemble contigs
    auto result = OLC(sequences, minOverlap);
    vector<string> contigs = result.first;
    vector<string> unassembledReads = result.second;

    // Print the contigs
    cout << "Contigs:" << endl;
    for (const auto& contig : contigs) {
        cout << contig << endl;
    }

    // Print unassembled reads
    cout << "Unassembled reads:" << endl;
    for (const auto& read : unassembledReads) {
        cout << read << endl;
    }

    return 0;
}