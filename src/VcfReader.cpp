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
        free(fwdCounts); // Free allocated memory
    }

    // Extract REV (reverse strand supporting reads)
    int32_t *revCounts = nullptr;
    int nrev = 0; // Initialize the number of elements to 0
    nrev = bcf_get_info_int32(header, record, "REV", &revCounts, &nrev);
    if (nrev > 0) {
        var.revCounts = *revCounts;
        free(revCounts); // Free allocated memory
    }

    return var;
}
