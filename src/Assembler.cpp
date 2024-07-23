#include <iostream>
#include <string> // carreguem la llibreria d'strings
#include <fstream> // llibreria per obrir inputs
#include <vector>
#include <list>
#include <algorithm>
#include <map>
#include "Assembler.h"
#include "ssw_cpp.h"
#include <unordered_map>
using namespace std;

// Implement the Assembler constructor
Assembler::Assembler(std::vector<std::string> reads, int kSize, int maxMismatches, bool debug)
    : reads(reads), kSize(kSize), maxMismatches(maxMismatches), debug(debug), totalAssembled(0), totalReads(reads.size()) {
    // Additional initialization if needed
}


// Function to calculate k-mer diversity
float kmerDiversity(const std::string& sequence, int k=3) {
    std::unordered_set<std::string> uniqueKmers;
    int n = sequence.size();

    // Extract all k-mers and add to the set
    for (int i = 0; i <= n - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        uniqueKmers.insert(kmer);
    }

    // Calculate diversity as the number of unique k-mers divided by the total possible k-mers
    int totalKmers = n - k + 1;
    if (totalKmers <= 0) return 0.0f; // Handle edge case where sequence length < k
    float diversity = static_cast<float>(uniqueKmers.size()) / totalKmers;

    return diversity;
}


string inline revComp(string seq) {

	string rev_seq(seq);
	reverse(rev_seq.begin(),rev_seq.end());
	for (int i = 0; i < rev_seq.length(); i++) {
		switch (rev_seq[i]) 
		{
			case 'A': rev_seq[i] = 'T'; break;
			case 'T': rev_seq[i] = 'A'; break;
			case 'C': rev_seq[i] = 'G'; break;
			case 'G': rev_seq[i] = 'C'; break;
			case 'a': rev_seq[i] = 't'; break;
			case 't': rev_seq[i] = 'a'; break;
			case 'c': rev_seq[i] = 'g'; break;
			case 'g': rev_seq[i] = 'c'; break;
		}
	}
	return rev_seq; 
}

std::vector<Nucleotide> revertConsensus( vector<Nucleotide>& vec) {

	std::vector<Nucleotide> outVec;
	for (auto& b : vec) {
		Nucleotide ntd;
		ntd.A = b.T;
		ntd.C = b.G;
		ntd.G = b.C;
		ntd.T = b.A;
		ntd.totalBases = b.totalBases;
		outVec.push_back(ntd);
	}
	reverse(outVec.begin(),outVec.end());
	return outVec;
}

std::multimap<std::string, std::string> Assembler::createReadTable() {

	std::multimap<std::string, std::string> readHash;

	for (auto read : reads) {
		readHash.insert(std::pair<std::string, std::string>(read, read));
	}
	return readHash;
}

// std::unordered_multimap<std::string, Contig> Assembler::createPrefixTable() {
//     std::unordered_multimap<std::string, Contig> prefixHash;

//     for (auto contig : ContigList) {
//         for (int k = 15; k <= contig.seq.length(); ++k) { // Use k-mers of increasing length from 15 to full length
//             std::string prefix = contig.seq.substr(0, k);
//             std::string revCompRead = revComp(contig.seq);
//             std::string revCompPrefix = revCompRead.substr(0, k);
//             Contig revCompContig = contig;
//             revCompContig.seq = revCompRead;

//             prefixHash.insert(std::pair<std::string, Contig>(prefix, contig));
//             prefixHash.insert(std::pair<std::string, Contig>(revCompPrefix, contig));
//         }
//     }

//     return prefixHash;
// }


std::unordered_multimap<std::string, Contig> Assembler::createPrefixTable() {

	std::unordered_multimap<std::string, Contig> prefixHash;

	for (auto contig : ContigList) {

		// Equal than SSAKE, VCAKE, popuating the first 11 bases on a hash table with the value being the full sequence
		std::string prefix = contig.seq.substr(0, kSize);
		std::string revCompRead = revComp(contig.seq);

		Contig revCompContig = contig;	 // Create revcomp contig struct
		revCompContig.seq = revCompRead; // RevComp of the contig
		std::string revCompPrefix = revCompRead.substr(0, kSize); // Prefix of revcomp

		prefixHash.insert(std::pair<std::string, Contig>(prefix, contig));
		prefixHash.insert(std::pair<std::string, Contig>(revCompPrefix, contig));
	}
	return prefixHash;
}

std::pair<std::string, int> mostVotedBase(Nucleotide Base) {
	int max = 0;
	std::string ntd;

	if (Base.A > max) {
		max = Base.A;
		ntd = "A";
	}
	if (Base.C > max) {
		max = Base.C;
		ntd = "C";
	}
	if (Base.T > max) {
		max = Base.T;
		ntd = "T";
	}
	if (Base.G > max) {
		max = Base.G;
		ntd = "G";
	}
	//cout << endl;
	return std::make_pair(ntd, max);	
}

int Assembler::getNumAssembled() {
	return totalAssembled;
}
int Assembler::getTotalReads() {
	return totalReads;
}

std::vector<std::string> Assembler::ungappedGreedy() {
    for (auto& r : reads) {
        Contig c;
        c.seq = r;
        for (auto& base : r) {
            Nucleotide ntd;
            ntd.totalBases++;
            switch (base) {
                case 'A': ntd.A = 1; break;
                case 'C': ntd.C = 1; break;
                case 'T': ntd.T = 1; break;
                case 'G': ntd.G = 1; break;
            }
            c.Consensus.push_back(ntd);
        }
        ContigList.push_back(c);
    }

    int notExtended = 0;
    std::map<std::string, int> seenMap;
    std::map<std::string, int> seenPair;
    std::multimap<std::string, std::string> readsHash = createReadTable();
    std::unordered_multimap<std::string, Contig> prefixTable = createPrefixTable();
	std::vector<std::string> testOut;
    // return testOut;

    int totalReads = readsHash.size();
    int round = 0;

    while (true) {
		// std::cout << "bucle" << std::endl;
        if (ContigList.empty()) break;

        Contig read = ContigList[0];
        round++;
        seenPair[read.seq]++;

        if (seenPair[read.seq] > 1) {
            bool allSeen = std::all_of(ContigList.begin(), ContigList.end(), [&seenPair](const Contig& c) {
                return seenPair[c.seq] > 1;
            });
            if (allSeen) break;
        }

        if (reads.size() <= 1) break;

        Contig contig;
        bool isRevComp = false;
        std::tie(notExtended, contig) = Extend(prefixTable, readsHash, seenPair, read, notExtended, isRevComp);

        isRevComp = true;
        std::tie(notExtended, contig) = Extend(prefixTable, readsHash, seenPair, contig, notExtended, isRevComp);
    }

    if (reads.size() > 1) {
        for (auto& r : readsHash) {
            reads.erase(std::remove(reads.begin(), reads.end(), r.first), reads.end());
        }
    }
    totalAssembled = totalReads - reads.size();
    std::vector<std::string> outReads;
    for (auto& r : reads) {
        if (r.length() >= kSize) {
            outReads.push_back(r);
        }
    }
    return outReads;
}




std::pair<int, Contig> Assembler::Extend( std::unordered_multimap<std::string, Contig>& prefixTable, std::multimap<std::string, std::string>& readHash, std::map<std::string, int>& seenPair, Contig& read, int notExtended, bool isRevComp ) {

	std::vector<candidateRead> arrayHits;
	std::vector<Nucleotide> consensusVector;
	std::string contig;
	bool isFound = false;
	std::string SEQ = read.seq;

	Contig outContig;

	if (SEQ.length() < kSize) {
		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
		reads.erase(std::remove(reads.begin(), reads.end(), SEQ), reads.end());
		reads.push_back(read.seq);
		ContigList.push_back(read);
		return std::make_pair(notExtended, outContig);
	}

	std::string revCompRead;
	if (isRevComp) {
		revCompRead = revComp(read.seq);
	}

	int min = 100000;
	int max = 0;
	int startOffset;
	int endOffset;
	int lastPos = read.seq.length();
	int span;
	if (isRevComp) {
		SEQ = revComp(read.seq);
	}

	int bestKSize = kSize;
	bool suitableKSizeFound = false;
	std::string testSeed = SEQ.substr(0, bestKSize);
	float seedComplexity = kmerDiversity(testSeed, 2);
	if (seedComplexity >= 0.1) {
		// bestKSize = kSize;
		suitableKSizeFound = true;
	}


	// for (int currentKSize = 15; currentKSize <= SEQ.length(); ++currentKSize) {
	// 	std::string testSeed = SEQ.substr(0, currentKSize);
	// 	float seedComplexity = kmerDiversity(testSeed, 2);
	// 	// std::cout << testSeed << " complexity:" << seedComplexity << std::endl;
	// 	if (seedComplexity >= 0.1) {
	// 		bestKSize = currentKSize;
	// 		suitableKSizeFound = true;
	// 		// if (debug) {
	// 		//     std::cout << "Optimized k-mer size: " << bestKSize << std::endl;
	// 		// }
	// 		break;
	// 	}
	// }

	if (!suitableKSizeFound) {
		// Remove the read from the assembly process
		// std::cout << "No suitable k-mer size found for read. Removing from assembly." << std::endl;
		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
		reads.erase(std::remove(reads.begin(), reads.end(), SEQ), reads.end());
		return std::make_pair(notExtended, outContig); // Return early
	}


// int bestKSize = 15;
// for (int currentKSize = 15; currentKSize <= SEQ.length(); ++currentKSize) {
//     std::string testSeed = SEQ.substr(0, currentKSize);
//     float seedComplexity = kmerDiversity(testSeed, 2);
//     std::cout << testSeed << " complexity:" << seedComplexity << std::endl;
//     if (seedComplexity >= 0.1) {
//         bestKSize = currentKSize;
//         if (debug) {
//             std::cout << "Optimized k-mer size: " << bestKSize << std::endl;
//         }
//         break;
//     }
    
//     // If we reach the end of the loop without finding a suitable kSize, use the last value
//     if (currentKSize == SEQ.length()) {
//         bestKSize = currentKSize;
//     }
// }
// kSize = bestKSize;
// std::cout << "Selected kSize: " << kSize << std::endl;

	for (int i = 0; i < SEQ.length()-kSize; i++) {

		std::string seed = SEQ.substr(i, kSize);
			
		std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret;
		ret = prefixTable.equal_range(seed);

		for (std::unordered_multimap<std::string, Contig>::iterator it=ret.first; it!=ret.second; ++it) {

			int matches = 0;
			int mismatches = 0;
			int m = 0;
			std::string SEQ2 = it->second.seq;

			if (isRevComp) {
				SEQ2 = revComp(SEQ2);
			}

   		    std::string SEQ_slice = SEQ.substr(i, SEQ.length());
       		std::string SEQ2_slice = SEQ2.substr(0, SEQ_slice.length());

			for (int n = 0; n < SEQ_slice.length(); n++) {
				if (SEQ_slice[n] == SEQ2_slice[n]) {
					matches++;
				}
				else {
					mismatches++;
				}
			}
			// int mmAllowd = (m * maxMismatches)/100;
			int mmAllowd = 1;
			if (mismatches <= mmAllowd ) {
				if (SEQ2 == SEQ) {					
					continue;
				}

				isFound = true;
				candidateRead cr;

				if (isRevComp) {
					cr.Seed = revComp(it->first);
					cr.Sequence = SEQ2;
					cr.ContigStruct   = it->second;
					cr.ContigStruct.seq = it->second.seq;
					std::vector<Nucleotide> tmpConsensus = it->second.Consensus;
					reverse(tmpConsensus.begin(),tmpConsensus.end());
					cr.Consensus = tmpConsensus;
				}
				else {				
					cr.Seed = it->first;
					cr.Sequence = SEQ2;
					cr.ContigStruct   = it->second;
					cr.Consensus =  it->second.Consensus;
				}
				cr.offset   = i;
				if (debug) {
					std::string spacer = " ";
					for (int l = 0; l < i; l++) {
						spacer += " ";
					}
					std::cout << " ---------------------" << std::endl;		
					std::cout <<  " " << SEQ << std::endl;
					std::cout << spacer << SEQ2 << std::endl;
				}

				// Debug: Print the slices of nucleotides to be compared
				// std::cout << "Comparing slices:\n";
				// std::cout << " mismatches:" << mismatches << std::endl;
				// std::cout << " allowed mismatches:" << mmAllowd << std::endl;
				// std::cout << " ---------------------" << std::endl;
				arrayHits.push_back(cr);
			}
		}
	}
	if (isFound==false) {

		notExtended++;
		contig = SEQ;

		if (isRevComp) {
			contig = revComp(contig);
			reads.erase(std::remove(reads.begin(), reads.end(), contig), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revComp(contig)), reads.end());

			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
			std::string prefixRead = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));
		}

		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
		reads.erase(std::remove(reads.begin(), reads.end(), contig), reads.end());

		// Do not remove the read itself if it cannot be extended, instead, append it at the end
		ContigList.push_back(read);
		reads.push_back(contig);
		std::string prefix = contig.substr(0, kSize);

		std::string revCompContig = revComp(contig);
		std::string prefRevComp = revCompContig.substr(0, kSize);

		outContig = read;
	}
	else {

		// Adding read offsets and others
		candidateRead cr;
		cr.Sequence  = SEQ;

		if (isRevComp) {
			std::vector<Nucleotide> tmpConsensus = read.Consensus;
			reverse(tmpConsensus.begin(),tmpConsensus.end());
			cr.Consensus = tmpConsensus;
		}
		else {
			cr.Consensus = read.Consensus;
		}
		cr.offset    = 0;
		cr.Seed      = read.seq.substr(0, kSize);
		arrayHits.push_back(cr);

		for (auto& h : arrayHits) {

			if (h.offset < min) {
				startOffset = h.offset;
				min = startOffset;
			}
			if (h.offset+ h.Sequence.length() > max) {
				endOffset = h.offset+ h.Sequence.length();
				max = endOffset;
			}
			std::string revCompRead = revComp(h.Sequence);

			if (h.Sequence != SEQ) {
				readHash.erase(h.Sequence);
				readHash.erase(revCompRead);
			}

			std::string prefixInitRead = h.Sequence.substr(0, kSize);
			std::string revCompPrefix = revCompRead.substr(0, kSize);


			std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret2;
			ret2 = prefixTable.equal_range(prefixInitRead);

			for (std::unordered_multimap<std::string, Contig>::iterator it=ret2.first; it!=ret2.second; ++it) {
				if (h.Sequence == it->second.seq or h.Sequence == revComp(it->second.seq)) {
					prefixTable.erase(it);
					break;
				}
			}

			std::pair <std::unordered_multimap<std::string, Contig>::iterator, std::unordered_multimap<std::string, Contig>::iterator> ret3;
			ret3 = prefixTable.equal_range(revCompPrefix);

			for (std::unordered_multimap<std::string, Contig>::iterator it=ret3.first; it!=ret3.second; ++it) {
				if (h.Sequence == it->second.seq or h.Sequence == revComp(it->second.seq)) {
					prefixTable.erase(it);		
					break;
				}
			}

			reads.erase(std::remove(reads.begin(), reads.end(), h.Sequence), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());

			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());
			ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), h.ContigStruct), ContigList.end());
		}

		// Initializing consesnsusVector
		for (int pos = 0; pos < endOffset; pos++) {
			Nucleotide base;
			consensusVector.push_back(base);
		}
		int numk = 0;
		int numHit = 0;
        int totalMismatches = 0;
        int mismatchThreshold = 5;  // Define a threshold for mismatches

		for(auto& r : arrayHits) {
			Nucleotide base;
			int j = r.offset;
			int i = 0;
			std::string spacer = "";
			if (debug) {
				for (int idx=0; idx<r.offset;idx++) {
					spacer+= " ";
				}
				std::cout << " " << " " << spacer << r.Sequence << std::endl;
			}

			std::string strout = "";
  			for (auto& b : r.Sequence){ 
				consensusVector[j].totalBases+= r.Consensus[i].totalBases;

				if (isRevComp) {
					consensusVector[j].A += r.Consensus[i].T;
					consensusVector[j].C += r.Consensus[i].G;
					consensusVector[j].T += r.Consensus[i].A;
					consensusVector[j].G += r.Consensus[i].C;
				}	
				else {
					consensusVector[j].A += r.Consensus[i].A;
					consensusVector[j].C += r.Consensus[i].C;
					consensusVector[j].T += r.Consensus[i].T;
					consensusVector[j].G += r.Consensus[i].G;
				}

                if (b != mostVotedBase(consensusVector[j]).first[0]) {
                    totalMismatches++;
                }

				switch(b) { 		
					case 'A':
						consensusVector[j].A++;
						break;
					case 'C':
						consensusVector[j].C++;
						break;		
					case 'T':
						consensusVector[j].T++;
						break;
					case 'G':
						consensusVector[j].G++;
						break;															
				}
				j++;
				i++;
			}
			numk = j;
			numHit++;
		}

        // if (totalMismatches > mismatchThreshold) {
        //     std::cout << "Contig discarded due to high mismatches: " << totalMismatches << std::endl;
      	// 	std::cout << std::endl;

        //     return std::make_pair(notExtended, outContig);
        // }

		// std::cout << std::endl;

		int x = 0;
		std::vector<Nucleotide> NewConsensus;

		for (auto& ntd : consensusVector) {
			std::string b;
			int totalb;
			if (ntd.totalBases > 0) {
				std::tie(b, totalb) = mostVotedBase(ntd);
				contig+=b;
				NewConsensus.push_back(ntd);
			}
		}
		if (isRevComp) {
			contig = revComp(contig);
		}

		ContigList.erase(std::remove(ContigList.begin(), ContigList.end(), read), ContigList.end());

		readHash.erase(read.seq);
		reads.erase(std::remove(reads.begin(), reads.end(), read.seq), reads.end());

		reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());
		readHash.erase(revCompRead);

		auto it = reads.insert(reads.begin(), contig);
		std::string prefixContig = contig.substr(0, kSize);

		Contig NewContig;
		NewContig.seq = contig;
		NewContig.Consensus = NewConsensus;

		// if (debug) {
		// 	cout << "\n" << "Contig => " << contig << endl;
		// }

		if (isRevComp) {
			NewContig.Consensus = revertConsensus(NewConsensus);
		}
		if (contig.length() >= kSize) {
			// if (debug) {
			// 	cout << "size is ok " << contig.length() << "\t" << kSize << "\t" << reads.size() << endl;
			// }
			ContigList.insert(ContigList.begin(), NewContig);

			//cout << endl;
			std::string prefixRead = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));

			prefixTable.insert(std::pair<std::string, Contig>(prefixContig, NewContig));
			prefixTable.insert(std::pair<std::string, Contig>(revComp(prefixContig), NewContig));

			outContig = NewContig;
		}
		else {
			// if (debug) {
			// 	cout << "size is notok " << contig.length() << endl;
			// }
			std::string prefixRead        = read.seq.substr(0, kSize);
			std::string refCompPrefixRead = revComp(prefixRead);

			prefixTable.erase(prefixRead);
			prefixTable.erase(revComp(prefixRead));

			reads.erase(std::remove(reads.begin(), reads.end(), read.seq), reads.end());
			reads.erase(std::remove(reads.begin(), reads.end(), revCompRead), reads.end());

			readHash.erase(revComp(contig));
			readHash.erase(contig);

		}

	}
	return std::make_pair(notExtended, outContig);
}