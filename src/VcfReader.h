#ifndef VCFREADER_H
#define VCFREADER_H

#include <string>
#include <vector>
#include <htslib/vcf.h>
#include "variant_t.h"

class VcfReader {
public:
    VcfReader(const std::string &vcfFilename);
    ~VcfReader();
    
    std::vector<variant_t> readVariants();

private:
    htsFile *vcfFile;
    bcf_hdr_t *header;
    bcf1_t *record;

    variant_t parseVariant(bcf1_t *record);
};

#endif // VCFREADER_H
