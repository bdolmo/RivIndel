#ifndef VCFWRITER_H
#define VCFWRITER_H

#include <fstream>
#include <string>
#include "variant_t.h"

class VcfWriter {
    public:
        VcfWriter(const std::string& vcfFilename, const std::string& bamFilename);
        ~VcfWriter();
        void writeHeader();
        void writeVariant(const variant_t& var);
    private:
        std::ofstream vcfFile;
        std::string sampleName;
        samFile* samFilePtr;
        bam_hdr_t* bamHeader;
};

#endif // VCFWRITER_H
