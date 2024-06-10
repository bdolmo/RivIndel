#include <iostream>
#include <string>
#include <vector>
#include "BamRecord.h"

#ifndef VARIANT_T_H
#define VARIANT_T_H

struct variant_t {
    std::string chr;
    int pos;
    std::string ref;
    std::string alt;
    std::string alleleDepth;
    float mapQual;
    int readSupport;
    int meanDepth;
    int totalDepth;   // Total number of reads covering the variant position
    float alleleFrequency;
    float kmerDiversity; // New member to store k-mer diversity
    float flankingKmerDiversityUpstream; // New member for upstream k-mer diversity
    float flankingKmerDiversityDownstream; // New member for downstream k-mer diversity
    float gcContent; // New member to store GC content

    int longestHomopolymerRun; // New member for the longest homopolymer run
    int numberOfHomopolymerRuns; // New member for the number of homopolymer run

    int longestDintdRun;
    int numberOfDintd;
    int plusStrand;
    int minusStrand;
    double strandBias;
    double meanBaseQuality;
};

#endif // VARIANT_T_H
