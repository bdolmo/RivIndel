#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <set>

using namespace std;

struct Overlap {
    int dest;   // Index of destination read
    int length; // Length of the overlap
};

// Finds maximum overlap with up to two mismatches allowed
int findMaxOverlapWithMismatches(const string& a, const string& b, int min_length) {
    int max_valid_overlap = 0;
    for (int len = min(a.size(), b.size()); len >= min_length; --len) {
        int mismatches = 0;
        for (int i = 0; i < len; ++i) {
            if (a[a.size() - len + i] != b[i]) {
                mismatches++;
                if (mismatches > 2) break;
            }
        }
        if (mismatches <= 2) {
            max_valid_overlap = len;
            break;
        }
    }
    return max_valid_overlap;
}

// Builds an overlap graph from the reads
map<int, vector<Overlap>> buildOverlapGraph(const vector<string>& reads, int min_length) {
    map<int, vector<Overlap>> graph;
    for (int i = 0; i < reads.size(); ++i) {
        for (int j = 0; j < reads.size(); ++j) {
            if (i != j) {
                int overlap = findMaxOverlapWithMismatches(reads[i], reads[j], min_length);
                if (overlap >= min_length) { // Ensure the overlap length is at least min_length
                    graph[i].push_back({j, overlap});
                }
            }
        }
    }
    
    // Debugging output to visualize the overlap graph
    cout << "Overlap Graph:" << endl;
    for (const auto& node : graph) {
        cout << "Read " << node.first << " overlaps with:" << endl;
        for (const auto& edge : node.second) {
            cout << "  -> Read " << edge.dest << " with overlap length " << edge.length << endl;
        }
    }
    
    return graph;
}

// Greedily finds the longest contig by following the longest overlaps
vector<int> findLongestContigPath(const vector<string>& reads, const map<int, vector<Overlap>>& graph, set<int>& visited, int start) {
    vector<int> path;
    int current = start;

    while (visited.find(current) == visited.end()) {
        visited.insert(current);
        path.push_back(current);

        if (graph.find(current) == graph.end() || graph.at(current).empty()) break;

        auto max_overlap = *max_element(graph.at(current).begin(), graph.at(current).end(),
            [](const Overlap& a, const Overlap& b) {
                return a.length < b.length;
            });

        if (visited.count(max_overlap.dest)) break;
        current = max_overlap.dest;
    }

    return path;
}

// Derives consensus sequence from the contig path
string deriveConsensus(const vector<string>& reads, const vector<int>& path, const map<int, vector<Overlap>>& graph) {
    if (path.empty()) return "";
    vector<map<char, int>> nucleotide_votes;
    int current_position = 0;

    for (size_t i = 0; i < path.size(); ++i) {
        const string& read = reads[path[i]];
        int overlap_length = 0;
        if (i > 0) {
            const auto& overlaps = graph.at(path[i - 1]);
            auto it = find_if(overlaps.begin(), overlaps.end(), [&](const Overlap& o) {
                return o.dest == path[i];
            });
            if (it != overlaps.end()) {
                overlap_length = it->length;
            }
        }

        for (int j = 0; j < read.size(); ++j) {
            if (current_position + j >= nucleotide_votes.size()) {
                nucleotide_votes.push_back(map<char, int>());
            }
            nucleotide_votes[current_position + j][read[j]]++;
        }

        current_position += (read.size() - overlap_length);
    }

    string consensus;
    for (const auto& position_votes : nucleotide_votes) {
        char most_voted = '-';
        int max_votes = 0;
        for (const auto& vote : position_votes) {
            if (vote.second > max_votes) {
                most_voted = vote.first;
                max_votes = vote.second;
            }
        }
        consensus += most_voted;
    }

    return consensus;
}


// Main assembly function using OLC-based approach
vector<string> greedyAssembler(vector<string>& reads, int min_length) {
    vector<string> contigs;
    set<int> visited;
    auto graph = buildOverlapGraph(reads, min_length);

    for (int i = 0; i < reads.size(); ++i) {
        if (visited.find(i) == visited.end()) {
            auto path = findLongestContigPath(reads, graph, visited, i);
            if (!path.empty()) {
                string consensus = deriveConsensus(reads, path, graph);
                contigs.push_back(consensus);
            }
        }
    }

    return contigs;
}

// int main() {
//     vector<string> reads = {
//         "AGCTAGCCTAGG",
//         "TAGGAGCTTTAC",
//         "TTACGGCAT",
//         "GGCATCATGCAAGTT",
//         "CATGCAAGTT"
//     };

//     int min_overlap_length = 25;
//     vector<string> contigs = greedyAssembler(reads, min_overlap_length);
//     cout << "Final Contigs: " << endl;
//     for (const string& contig : contigs) {
//         cout << contig << endl;
//     }

//     return 0;
// }
