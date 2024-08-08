#include "VcfReader.h"
#include "VcfWriter.h"
#include "BamReader.h"
#include <omp.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <boost/math/distributions/hypergeometric.hpp>

using namespace boost::math;

// Function to perform Fisher's exact test
double FisherExactTest(int a, int b, int c, int d) {

    // Total number of items
    int n = a + b + c + d;

    // Number of items in each row and column
    int row1 = a + b;
    int row2 = c + d;
    int col1 = a + c;
    int col2 = b + d;

    // Create the hypergeometric distribution
    hypergeometric_distribution<> hg1(col1, row1, n);
    hypergeometric_distribution<> hg2(col2, row1, n);

    // Calculate the p-value (two-tailed)
    double p_value = 0.0;

    // Iterate over all possible values
    for (int x = std::max(0, row1 + col1 - n); x <= std::min(row1, col1); ++x) {
        double p = pdf(hg1, x) * pdf(hg2, row1 - x);
        if (p <= pdf(hg1, a) * pdf(hg2, b)) {
            p_value += p;
        }
    }

    return p_value;
}



// Function to compare variants and classify them as germline or somatic
void classififyVariants(const std::string& tumorVcfFile, const std::string& tumorBamFile, const std::string& normalVcfFile, const std::string& outputVcfFile) {
    VcfReader tumorReader(tumorVcfFile);
    VcfReader normalReader(normalVcfFile);

    std::vector<variant_t> tumorVariants = tumorReader.readVariants();
    std::vector<variant_t> normalVariants = normalReader.readVariants();

    std::unordered_map<std::string, variant_t> normalVariantsMap;
    for (const auto& var : normalVariants) {
        std::string key = var.chr + ":" + std::to_string(var.pos) + ":" + var.ref + ":" + var.alt;
        normalVariantsMap[key] = var;
    }

    VcfWriter vcfWriter(outputVcfFile, tumorBamFile); // Assuming tumor BAM file for header info
    vcfWriter.writeHeader();

    for (auto& tumorVar : tumorVariants) {
        std::string key = tumorVar.chr + ":" + std::to_string(tumorVar.pos) + ":" + tumorVar.ref + ":" + tumorVar.alt;
        bool isSomatic = false;
        double p_value = 1.0;

        if (normalVariantsMap.find(key) != normalVariantsMap.end()) {
            variant_t normalVar = normalVariantsMap[key];
            int tumorTotal = tumorVar.totalDepth;
            int normalTotal = normalVar.totalDepth;
            int tumorAlt = tumorVar.readSupport;
            int normalAlt = normalVar.readSupport;

            // Adjust counts to avoid division by zero
            tumorAlt = std::max(1, tumorAlt);
            tumorTotal = std::max(tumorAlt + 1, tumorTotal);
            normalAlt = std::max(1, normalAlt);
            normalTotal = std::max(normalAlt + 1, normalTotal);
            // std::cout <<key<< " " << tumorAlt << " " <<  tumorTotal- tumorAlt << " "  << normalAlt << " " <<  normalTotal - normalAlt << "\n";
            p_value = FisherExactTest(tumorAlt, tumorTotal - tumorAlt, normalAlt, normalTotal - normalAlt);
            // std::cout << p_value << std::endl;

            if (p_value < 0.05) {
                isSomatic = true;
            } else {
                isSomatic = false;
            }
        } else {
            isSomatic = true;
        }

        // Add the classification and p-value to the INFO field
        tumorVar.status = isSomatic ? "SOMATIC" : "GERMLINE";
        // tumorVar.pValue = p_value;

        vcfWriter.writeVariant(tumorVar);
    }
}