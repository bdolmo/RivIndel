#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <utility>
#include <fstream>
#include "ssw_cpp.h"

using namespace std;

struct Nucleotide {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
    int totalBases = 0;
};

struct Contig {
    std::string seq;
    std::vector<Nucleotide> Consensus;
    bool operator==(const Contig& other) const {
        return seq == other.seq;
    }
};

struct candidateRead {
    std::string Seed;
    std::string Sequence;
    Contig ContigStruct;
    std::vector<Nucleotide> Consensus;
    int offset;
};


class Assembler {
	public:
		std::vector<std::string> reads;
		std::vector<Contig> ContigList;
		int kSize;
		int totalAssembled;
		int totalReads;
		bool debug;
		int maxMismatches;

		Assembler(std::vector<std::string> reads, int kSize, int maxMismatches, bool debug = false);

		std::multimap<std::string, std::string> createReadTable();
		std::unordered_multimap<std::string, Contig> createPrefixTable();
		int getNumAssembled();
		int getTotalReads();
		std::vector<std::string> ungappedGreedy();
		std::pair<int, Contig> Extend(std::unordered_multimap<std::string, Contig>& prefixTable, 
									std::multimap<std::string, std::string>& readHash, 
									std::map<std::string, int>& seenPair, 
									Contig& read, int notExtended, bool isRevComp);
};

// Function prototypes
string revComp(string seq);
std::vector<Nucleotide> revertConsensus(std::vector<Nucleotide>& vec);
std::pair<std::string, int> mostVotedBase(Nucleotide Base);

#endif // ASSEMBLER_HPP
