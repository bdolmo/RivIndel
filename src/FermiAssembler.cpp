#include "fml.h"
#include <vector>
#include <string>
#include <cstring>

// New function to perform assembly using fermi-lite with error correction
std::vector<std::string> assembleWithFermiLite(const std::vector<std::string>& reads) {
    std::vector<std::string> contigs;

    int n_seqs = reads.size();
    bseq1_t *seqs = new bseq1_t[n_seqs]();
    for (int i = 0; i < n_seqs; ++i) {
        seqs[i].l_seq = reads[i].length();
        seqs[i].seq = strdup(reads[i].c_str());
        
        seqs[i].qual = new char[reads[i].length() + 1];
        memset(seqs[i].qual, 'I', reads[i].length());
        seqs[i].qual[reads[i].length()] = '\0';
    }

    // Initialize default options
    fml_opt_t opt;

    fml_opt_init(&opt);

    // Set custom parameters
    // opt.ec_k = 11;        // Set k-mer length for error correction
    // opt.min_asm_ovlp = 11; // Set minimum overlap length during assembly
    // opt.min_merge_len = 11; // Set minimum length for merging overlaps

    // Perform error correction
    float kcov = fml_correct(&opt, n_seqs, seqs);
    // Optional: adjust parameters based on corrected sequences
    fml_opt_adjust(&opt, n_seqs, seqs);



    // Perform assembly
    int n_utgs;
    fml_utg_t *utgs = fml_assemble(&opt, n_seqs, seqs, &n_utgs);

    for (int i = 0; i < n_utgs; ++i) {
        std::string contig = std::string(utgs[i].seq, utgs[i].len);
        contigs.push_back(contig);
    }

    // It is not necessary to free "seqs" allocated memory, since fermi-lite handles its ownership and removal internally!!
    fml_utg_destroy(n_utgs, utgs);

    return contigs;
}

// #include "fml.h"

// // New function to perform assembly using fermi-lite
// std::vector<std::string> assembleWithFermiLite(const std::vector<std::string>& reads) {
//     std::vector<std::string> contigs;

//     int n_seqs = reads.size();
//     bseq1_t *seqs = new bseq1_t[n_seqs]();
//     for (int i = 0; i < n_seqs; ++i) {
//         seqs[i].l_seq = reads[i].length();
//         seqs[i].seq = strdup(reads[i].c_str());
        
//         seqs[i].qual = new char[reads[i].length() + 1];
//         memset(seqs[i].qual, 'I', reads[i].length());
//         seqs[i].qual[reads[i].length()] = '\0';
//     }
//     // Perform assembly
//     fml_opt_t opt;
//     fml_opt_init(&opt);
//     int n_utgs;

//     fml_utg_t *utgs = fml_assemble(&opt, n_seqs, seqs, &n_utgs);

//     for (int i = 0; i < n_utgs; ++i) {
//         std::string contig = std::string(utgs[i].seq, utgs[i].len);
//         contigs.push_back(std::string(utgs[i].seq, utgs[i].len));
//     }

//     // It is not necessary to free "seqs" allocated memory, since fermi-lite handles its ownership and removal internally!!
//     fml_utg_destroy(n_utgs, utgs);

//     return contigs;
// }