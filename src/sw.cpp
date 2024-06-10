#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>
#include <regex>
#include <set>
// Scoring functions with parameters
int score_match(char a, char b, int match, int mismatch) {
    return a == b ? match : mismatch;
}

int score_match_special(char a, char b, int match_special, int mismatch_special) {
    return a == b ? match_special : mismatch_special;
}

std::string compact_cigar_string(const std::string& cigar_str) {
    if (cigar_str.empty()) return ""; // Return an empty string if the input is empty

    std::string compact_cigar;
    int count = 1; // Start counting from 1 since we always have at least one operation

    for (size_t i = 1; i <= cigar_str.size(); ++i) { // Start from the second character
        if (i < cigar_str.size() && cigar_str[i] == cigar_str[i - 1]) {
            count++; // Increment count if the current and previous characters are the same
        } else {
            // Once we reach a different character or the end of the string,
            // append the count and the operation type to the compact string
            compact_cigar += std::to_string(count) + cigar_str[i - 1];
            count = 1; // Reset count for the next operation
        }
    }
    return compact_cigar;
}

struct AlignmentResult {
    int max_i;
    int max_j;
    int max_score;
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
    std::string seq1_align;
    std::string seq2_align;
    std::string spacer;
    std::string cigar;
    std::string cigarExtended;
    int non_match_operations;
};

AlignmentResult affine_semiglobal(
    const std::string& seq1, 
    const std::string& seq2, 
    int match, 
    int mismatch, 
    int gap_open, 
    int gap_extend, 
    int match_special, 
    int mismatch_special) 
{
    int m = seq1.length();
    int n = seq2.length();

    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> X(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> Y(m + 1, std::vector<int>(n + 1, 0));

    for (int i = 1; i <= m; ++i) {
        X[i][0] = 0;
        Y[i][0] = INT_MIN / 2; 
    }
    for (int j = 1; j <= n; ++j) {
        Y[0][j] = 0;
        X[0][j] = INT_MIN / 2;
    }

    // Filling matrices
    int max_score = 0, max_i = 0, max_j = 0;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            M[i][j] = std::max(
                0, 
                std::max({
                    M[i - 1][j - 1] + score_match(seq1[i - 1], seq2[j - 1], match, mismatch), 
                    X[i - 1][j - 1] + score_match_special(seq1[i - 1], seq2[j - 1], match_special, mismatch_special), 
                    Y[i - 1][j - 1] + score_match_special(seq1[i - 1], seq2[j - 1], match_special, mismatch_special),
                })
            );

            X[i][j] = std::max({
                0, 
                M[i - 1][j] + gap_open,
                X[i - 1][j] + gap_extend, 
                Y[i - 1][j]  
            });

            // Update rule for Y (deletions in seq1 or insertions in seq2)
            Y[i][j] = std::max({
                0, 
                M[i][j - 1] + gap_open,
                Y[i][j - 1] + gap_extend,
                X[i][j - 1]
            });

            // Update max score
            int current_max = std::max({M[i][j], X[i][j], Y[i][j]});
            if (current_max > max_score) {
                max_score = current_max;
                max_i = i;
                max_j = j;
            }
        }
    }
    int query_start = max_i;
    int query_end = max_i;
    int ref_start = max_j;
    int ref_end = max_j;
    int num_mismatches = 0;

    // Traceback
    std::string seq1_align, seq2_align, spacer, cigar_str;
    int i = max_i, j = max_j;
    while (i > 0 || j > 0) {
        int current_max = std::max({M[i][j], X[i][j], Y[i][j]});
        if (current_max == 0) break;

        if (current_max == M[i][j]) {
            seq1_align = seq1[i - 1] + seq1_align;
            seq2_align = seq2[j - 1] + seq2_align;
            spacer = (seq1[i - 1] == seq2[j - 1] ? '|' : '*') + spacer;

            if (seq1[i-1] == seq2[j-1]) {
                cigar_str = "M" + cigar_str;
            }
            else {
                cigar_str = "X" + cigar_str;
                num_mismatches++;
            }
            i--;
            j--;
        } else if (current_max == X[i][j]) {
            seq1_align = seq1[i - 1] + seq1_align;
            seq2_align = '-' + seq2_align;
            spacer = ' ' + spacer;
            cigar_str = "D" + cigar_str;
            i--;
        } else { 
            seq1_align = '-' + seq1_align;
            seq2_align = seq2[j - 1] + seq2_align;
            spacer = ' ' + spacer;
            cigar_str = "I" + cigar_str;
            j--;
        }
        if (i > 0) query_start = i - 1;
        if (j > 0) ref_start = j - 1;
    }


    std::string compact_cigar = compact_cigar_string(cigar_str);

    // Count the number of different operations apart from M and X
    std::vector<char> operations;
    for (char c : compact_cigar) {
        if (std::isalpha(c) && c != 'M') {
            operations.push_back(c);
        }
        
    }
    int non_match_operations = operations.size();

    return AlignmentResult{
        max_i, max_j, max_score,
        query_start, query_end,
        ref_start, ref_end,
        seq1_align, seq2_align, spacer,
        compact_cigar,
        cigar_str,
        non_match_operations
    };
}

// int main() {
//     std::string seq1 = "CCATCTCACAATTGCCAGTTAACGTCTTCCTTCTCTCTCTGTCATAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATTCCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTCGATGTGAGTTTCTGCTTTGCTGTGTGGGGGTCCATGGCTCTGAACCTCAGGCCCAC";
//     std::string seq2 = "GCTATCAAGGAATTAAGAGAAACCAACATCGATGTGAGTTTCTGCTTTGCTGTGTGGGGG";
    
//     // Parameters for scoring
//     int match = 2;
//     int mismatch = -6;
//     int gap_open = -18;
//     int gap_extend = 0;
//     int match_special = -12;
//     int mismatch_special = -20;

//     AlignmentResult result = affine_semiglobal(seq1, seq2, match, mismatch, gap_open, gap_extend, match_special, mismatch_special);
    
//     std::cout << "Alignment score: " << result.max_score << std::endl;
//     std::cout << "       " << result.seq1_align << std::endl;
//     std::cout << "       " << result.spacer << std::endl;
//     std::cout << "       " << result.seq2_align << std::endl;
//     std::cout << "Cigar: " << result.cigar << std::endl;

//     return 0;
// }