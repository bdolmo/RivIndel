#include "VcfReader.h"
#include <stdexcept>
#include <iostream>

VcfReader::VcfReader(const std::string &vcfFilename) : vcfFile(nullptr), header(nullptr), record(nullptr) {
    vcfFile = bcf_open(vcfFilename.c_str(), "r");
    if (!vcfFile) {
        throw std::runtime_error("Could not open VCF file: " + vcfFilename);
    }

    header = bcf_hdr_read(vcfFile);
    if (!header) {
        bcf_close(vcfFile);
        throw std::runtime_error("Could not read VCF header: " + vcfFilename);
    }

    record = bcf_init();
    if (!record) {
        bcf_hdr_destroy(header);
        bcf_close(vcfFile);
        throw std::runtime_error("Could not allocate VCF record");
    }
}

VcfReader::~VcfReader() {
    if (record) {
        bcf_destroy(record);
    }
    if (header) {
        bcf_hdr_destroy(header);
    }
    if (vcfFile) {
        bcf_close(vcfFile);
    }
}

std::vector<variant_t> VcfReader::readVariants() {
    std::vector<variant_t> variants;

    while (bcf_read(vcfFile, header, record) == 0) {
        bcf_unpack(record, BCF_UN_ALL); // Unpack all fields

        variant_t var = parseVariant(record);
        variants.push_back(var);
    }

    return variants;
}

variant_t VcfReader::parseVariant(bcf1_t *record) {
    variant_t var;
    var.chr = bcf_hdr_id2name(header, record->rid);
    var.pos = record->pos + 1; // VCF is 1-based

    // Extract REF and ALT alleles
    var.ref = std::string(record->d.allele[0]);
    if (record->n_allele > 1) {
        var.alt = std::string(record->d.allele[1]);
    } else {
        var.alt = ".";
    }

    // Extract additional INFO fields as needed
    // Example for DP (depth)
    int32_t *dp = nullptr;
    int ndp = 0; // Initialize the number of elements to 0
    ndp = bcf_get_info_int32(header, record, "DP", &dp, &ndp);
    if (ndp > 0) {
        var.totalDepth = *dp;
        free(dp); // Free allocated memory
    }

    // Extract AF (allele frequency)
    int32_t *ac = nullptr;
    int nac = 0; // Initialize the number of elements to 0
    nac = bcf_get_info_int32(header, record, "AC", &ac, &nac);
    if (nac > 0) {
        var.readSupport = *ac;
        free(ac); // Free allocated memory
    }


    // Extract AF (allele frequency)
    float *af = nullptr;
    int naf = 0; // Initialize the number of elements to 0
    naf = bcf_get_info_float(header, record, "AF", &af, &naf);
    if (naf > 0) {
        var.alleleFrequency = *af;
        free(af); // Free allocated memory
    }

    // Extract FWD (forward strand supporting reads)
    int32_t *fwdCounts = nullptr;
    int nfwd = 0; // Initialize the number of elements to 0
    nfwd = bcf_get_info_int32(header, record, "FWD", &fwdCounts, &nfwd);
    if (nfwd > 0) {
        var.fwdCounts = *fwdCounts;
        var.plusStrand = *fwdCounts;
        free(fwdCounts); // Free allocated memory
    }

    // Extract REV (reverse strand supporting reads)
    int32_t *revCounts = nullptr;
    int nrev = 0; // Initialize the number of elements to 0
    nrev = bcf_get_info_int32(header, record, "REV", &revCounts, &nrev);
    if (nrev > 0 && revCounts != nullptr) {
        var.revCounts = *revCounts;
        var.minusStrand = *revCounts;
        free(revCounts); // Free allocated memory
    }


    float *mapQual = nullptr; 
    int nmapQual = 0;  // Initialize the number of elements to 0
    // Retrieve the KDIV_CONTIG info as floats
    nmapQual = bcf_get_info_float(header, record, "MQ", &mapQual, &nmapQual);
    if (nmapQual > 0 && mapQual != nullptr) {
        // If nkmerdiv > 0, store the first value (or handle as needed)
        var.mapQual = *mapQual;
        free(mapQual); // Free allocated memory
    }



    float *kmerDivContig = nullptr;  // Pointer to store k-mer diversity
    int nkmerdiv = 0;  // Initialize the number of elements to 0
    // Retrieve the KDIV_CONTIG info as floats
    nkmerdiv = bcf_get_info_float(header, record, "KDIV_CONTIG", &kmerDivContig, &nkmerdiv);
    if (nkmerdiv > 0 && kmerDivContig != nullptr) {
        // If nkmerdiv > 0, store the first value (or handle as needed)
        var.kmerDiversity = kmerDivContig[0];
        free(kmerDivContig); // Free allocated memory
    }


    float *kmerdivFlank5 = nullptr; // New member to store k-mer diversity
    int nkmerdivFlank5 = 0; // Initialize the number of elements to 0
    nkmerdivFlank5 = bcf_get_info_float(header, record, "KDIV_5Flank", &kmerdivFlank5, &nkmerdivFlank5);
    if (nkmerdivFlank5 > 0) {
        var.flankingKmerDiversityUpstream = *kmerdivFlank5;
        free(kmerdivFlank5); // Free allocated memory
    }


    float *kmerdivFlank3 = nullptr; // New member to store k-mer diversity
    int nkmerdivFlank3 = 0; // Initialize the number of elements to 0
    nkmerdivFlank3 = bcf_get_info_float(header, record, "KDIV_3Flank", &kmerdivFlank3, &nkmerdivFlank3);
    if (nkmerdivFlank3 > 0) {
        var.flankingKmerDiversityDownstream = *kmerdivFlank3;
        free(kmerdivFlank3); // Free allocated memory
    }


    float *gcContent = nullptr; // New member to store k-mer diversity
    int ngcContent = 0; // Initialize the number of elements to 0
    ngcContent = bcf_get_info_float(header, record, "GC", &gcContent, &ngcContent);
    if (ngcContent > 0) {
        var.gcContent = *gcContent;
        free(gcContent); // Free allocated memory
    }

    float *errorRate = nullptr;  // Allocate memory for float pointer
    int nerrorRate = 0;  // Initialize the number of elements to 0
    nerrorRate = bcf_get_info_float(header, record, "ERR", &errorRate, &nerrorRate);

    if (nerrorRate > 0 && errorRate != nullptr) {
        double err = 0.0;
        // Convert from float to double by summing and averaging the float array
        for (int i = 0; i < nerrorRate; ++i) {
            err += static_cast<double>(errorRate[i]);
        }
        err /= nerrorRate;  // Compute the average if needed
        var.errorRate = err;  // Store the result as a double
        free(errorRate);  // Free allocated memory
    }

    // Chimeric Rate
    float *chimericRate = nullptr;  // Allocate memory for float pointer
    int nchimericRate = 0;  // Initialize the number of elements to 0
    nchimericRate = bcf_get_info_float(header, record, "CHIMR", &chimericRate, &nchimericRate);

    if (nchimericRate > 0 && chimericRate != nullptr) {
        double chimRate = 0.0;

        // Convert from float to double by summing and averaging the float array
        for (int i = 0; i < nchimericRate; ++i) {
            chimRate += static_cast<double>(chimericRate[i]);
        }
        chimRate /= nchimericRate;  // Compute the average if needed

        var.chimericRate = chimRate;  // Store the result as a double

        free(chimericRate);  // Free allocated memory
    }

    // Soft Clipped Rate
    float *softClippedRate = nullptr;  // Allocate memory for float pointer
    int nsoftClippedRate = 0;  // Initialize the number of elements to 0
    nsoftClippedRate = bcf_get_info_float(header, record, "SCR", &softClippedRate, &nsoftClippedRate);

    if (nsoftClippedRate > 0 && softClippedRate != nullptr) {
        double scrRate = 0.0;

        // Convert from float to double by summing and averaging the float array
        for (int i = 0; i < nsoftClippedRate; ++i) {
            scrRate += static_cast<double>(softClippedRate[i]);
        }
        scrRate /= nsoftClippedRate;  // Compute the average if needed

        var.softClippedRate = scrRate;  // Store the result as a double

        free(softClippedRate);  // Free allocated memory
    }

    // Strand Bias
    float *strandBias = nullptr;  // Allocate memory for float pointer
    int nstrandBias = 0;  // Initialize the number of elements to 0
    nstrandBias = bcf_get_info_float(header, record, "SB", &strandBias, &nstrandBias);

    if (nstrandBias > 0 && strandBias != nullptr) {
        double sb = 0.0;

        // Convert from float to double by summing and averaging the float array
        for (int i = 0; i < nstrandBias; ++i) {
            sb += static_cast<double>(strandBias[i]);
        }
        sb /= nstrandBias;  // Compute the average if needed

        var.strandBias = sb;  // Store the result as a double

        free(strandBias);  // Free allocated memory
    }



    float *meanBaseQuality = nullptr;  // New member to store the float values
    int nmeanBaseQuality = 0;  // Initialize the number of elements to 0
    // Retrieve the "MBQ" information as floats
    nmeanBaseQuality = bcf_get_info_float(header, record, "MBQ", &meanBaseQuality, &nmeanBaseQuality);
    if (nmeanBaseQuality > 0 && meanBaseQuality != nullptr) {
        double mbq = 0.0;
        // Convert from float to double by summing and averaging the float array
        for (int i = 0; i < nmeanBaseQuality; ++i) {
            mbq += static_cast<double>(meanBaseQuality[i]);
        }
        mbq /= nmeanBaseQuality;  // Compute the average if needed
        var.meanBaseQuality = mbq;  // Store the result as a double
        free(meanBaseQuality);  // Free allocated memory
    }



    // Extract REV (reverse strand supporting reads)
    int32_t *longestHomopolymerRun = nullptr;
    int nlongestHomopolymerRun = 0; // Initialize the number of elements to 0
    nlongestHomopolymerRun = bcf_get_info_int32(header, record, "HOML", &longestHomopolymerRun, &nlongestHomopolymerRun);
    if (nlongestHomopolymerRun > 0) {
        var.longestHomopolymerRun = *longestHomopolymerRun;
        free(longestHomopolymerRun); // Free allocated memory
    }

    // Extract REV (reverse strand supporting reads)
    int32_t *numberOfHomopolymerRuns = nullptr;
    int nnumberOfHomopolymerRuns = 0; // Initialize the number of elements to 0
    nnumberOfHomopolymerRuns = bcf_get_info_int32(header, record, "HOMN", &numberOfHomopolymerRuns, &nnumberOfHomopolymerRuns);
    if (nnumberOfHomopolymerRuns > 0) {
        var.numberOfHomopolymerRuns = *numberOfHomopolymerRuns;
        free(numberOfHomopolymerRuns); // Free allocated memory
    }

    // Extract REV (reverse strand supporting reads)
    int32_t *longestDintdRun = nullptr;
    int nlongestDintdRun = 0; // Initialize the number of elements to 0
    nlongestDintdRun = bcf_get_info_int32(header, record, "DNTDL", &longestDintdRun, &nlongestDintdRun);
    if (nlongestDintdRun > 0) {
        var.longestDintdRun = *longestDintdRun;
        free(longestDintdRun); // Free allocated memory
    }


    // Extract REV (reverse strand supporting reads)
    int32_t *numberOfDintd = nullptr;
    int nnumberOfDintd = 0; // Initialize the number of elements to 0
    nnumberOfDintd = bcf_get_info_int32(header, record, "DNTDN", &numberOfDintd, &nnumberOfDintd);
    if (nnumberOfDintd > 0) {
        var.numberOfDintd = *numberOfDintd;
        free(numberOfDintd); // Free allocated memory
    }


    // int plusStrand;
    // int minusStrand;


    return var;
}
