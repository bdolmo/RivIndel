#include <iostream>
#include <string>
#include <vector>
#include "BamRecord.h"

struct clustered_aln_t {
    std::string chromosome;
    int pos;
    BamRecord read;
    std::string bucket;
    char strand;
    double mean_bq = 0.0;
};