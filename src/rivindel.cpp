#include <iostream>
#include <string>
#include <htslib/sam.h>
#include "BamRecord.h"
#include "BamReader.h"
#include <map>
#include <vector>
#include <optional>
#include "Cluster.cpp"
#include "Detect.cpp"
#include "VcfWriter.h"
#include "fml.h"
#include "variant_t.h"


// #include "cluster_t.h"

int main(int argc, char *argv[]) {

    std::string bamFile;
    std::string refFasta;
    std::string vcfOut;
    std::string targetsBed;

    int minOpSize = 1;
    int maxSeparation = 5; 
    int minSupport = 3;
    int numThreads = 1;


    // Basic command-line argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--bam" && i + 1 < argc) {
            bamFile = argv[++i];
        } 
        else if (arg == "--ref" && i + 1 < argc) {
            refFasta = argv[++i];
        }
        else if (arg == "--bed" && i + 1 < argc) {
            targetsBed = argv[++i];
        }
        else if (arg == "--vcf" && i + 1 < argc) {
            vcfOut = argv[++i];
        } 
        else if (arg == "--minOpSize" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> minOpSize;
        } 
        else if (arg == "--maxSeparation" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> maxSeparation;
        } 
        else if (arg == "--minSupport" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> minSupport;
        } 
        else if (arg == "--threads" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> numThreads;
        } 

        else {
            std::cerr << " ERROR: Unknown argument or missing value for: " << arg << std::endl;
            return 1;
        }
    }

    if (bamFile.empty() || refFasta.empty() || vcfOut.empty()) {
        std::cerr << "Usage: " << "./rivindel" << std::endl << "Version:." << std::endl << std::endl <<
                "Input parameters" << std::endl <<
                "   --bam   STRING  Input BAM file" << std::endl <<
                "   --ref   STRING  Reference genome in FASTA format" << std::endl <<
                "   --vcf   STRING  Output VCF file" << std::endl <<
                "   --threads   INTEGER Number of threads (default=1)" << std::endl << std::endl <<

                "For targeted sequencing:" << std::endl <<
                "   --bed   STRING  Targeted regions BED file" << std::endl << std::endl <<

                "Thresholds: " << std::endl <<
                "   --minOpSize INTEGER Minimum number of CIGAR operations to consider a complex Indel" << std::endl <<
                "   --maxSeparation INTEGER Minimum separation in bp between CIGAR operations" << std::endl <<
                "   --minSupport    INTEGER Minimum read support" << std::endl;
        return 1;
    }

    std::filesystem::path bamPath(bamFile);
    std::string baseName = bamPath.stem().string(); // Get the basename without extension
    std::string tempBamFile = bamPath.parent_path().string() + "/" + baseName + ".rivindel.signals.bam"; // Appe

    // identify cluster reads with compatible signals for cxindels
    std::map<std::string, std::vector<clustered_aln_t>> clusteredReads = clusterInformativeReads(bamFile, 
        minOpSize, maxSeparation, minSupport, tempBamFile, targetsBed);


    // Assemble clusters, align to reference and call variants
    std::map<std::string, variant_t> complexVariants = DetectComplexIndels(clusteredReads, 
        bamFile, refFasta, numThreads);

    PopulateVCF(complexVariants, bamFile, vcfOut);

    return 0;
}