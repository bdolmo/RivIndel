#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>

// Helper function to calculate the reverse complement of a DNA sequence
std::string reverse_complement(const std::string& seq) {
    std::string rev_seq(seq.rbegin(), seq.rend());
    std::transform(rev_seq.begin(), rev_seq.end(), rev_seq.begin(), [](char c) {
        switch (c) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
        }
        return 'N';
    });
    return rev_seq;
}

// Function to populate bin and set hash tables
void populate_hash_tables(const std::vector<std::string>& reads,
                          std::unordered_map<std::string, int>& bin,
                          std::unordered_map<std::string, std::vector<std::string>>& set) {
    for (const auto& read : reads) {
        std::string key = read.substr(0, 11);  // First 11 bases as key
        std::string rev_read = reverse_complement(read);
        std::string rev_key = rev_read.substr(0, 11);

        bin[read] += 1;
        bin[rev_read] += 1;

        set[key].push_back(read);
        set[rev_key].push_back(rev_read);
    }
}

// Function to extend a seed sequence into a contig in one direction
std::string extend_contig(const std::string& seed, const std::unordered_map<std::string, int>& bin,
                          const std::unordered_map<std::string, std::vector<std::string>>& set, 
                          int t, int c) {
    std::string contig = seed;
    std::string current_suffix = contig.substr(contig.size() - 11, 11);

    while (true) {
        std::unordered_map<char, int> base_votes;
        int max_vote = 0;
        char best_base = 0;

        if (set.count(current_suffix)) {
            for (const std::string& read : set.at(current_suffix)) {
                if (read.size() > 11) {
                    char candidate_base = read[11];
                    base_votes[candidate_base] += bin.at(read);
                    if (base_votes[candidate_base] > max_vote) {
                        max_vote = base_votes[candidate_base];
                        best_base = candidate_base;
                    }
                }
            }
        }

        if (max_vote >= t && best_base != 0) {
            contig += best_base;
            current_suffix = contig.substr(contig.size() - 11, 11);
        } else {
            break;
        }
    }

    return contig;
}

// Main VCAKE algorithm function
std::vector<std::string> vcake(const std::vector<std::string>& reads, int n, int t, int m, int e, int c, int v, int x) {
    std::unordered_map<std::string, int> bin;
    std::unordered_map<std::string, std::vector<std::string>> set;
    std::vector<std::string> contigs;

    populate_hash_tables(reads, bin, set);

    for (const auto& entry : set) {
        std::string seed = entry.second.front();  // Use the first read as seed
        std::string forward_contig = extend_contig(seed, bin, set, t, c);
        std::string reverse_seed = reverse_complement(seed);
        std::string reverse_contig = extend_contig(reverse_seed, bin, set, t, c);

        std::string full_contig = reverse_complement(reverse_contig) + forward_contig.substr(11);
        if (!full_contig.empty()) {
            contigs.push_back(full_contig);
        }
    }

    return contigs;
}

// int main() {
//     std::vector<std::string> reads = {"ACGTACGTACGT", "CGTACGTACGTA", "GTACGTACGTAC", "TACGTACGTACG"};
//     int n = 15, t = 1, m = 12, e = 10, c = 1, v = 2, x = 10;
//     std::vector<std::string> result = vcake(reads, n, t, m, e, c, v, x);
//     for (const std::string& contig : result) {
//         std::cout << "Contig: " << contig << std::endl;
//     }
//     return 0;
// }
