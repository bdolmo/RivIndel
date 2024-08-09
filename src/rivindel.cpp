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
#include <omp.h>
#include <set>
#include <filesystem>
#include <mutex>
#include "PairedAnalysis.cpp"

// Function declarations
void extractSignals(const std::string &bamFile, const std::vector<std::string> &chromosomes, const std::string &targetsBed, const std::string &excludeBed, int minOpSize, int maxSeparation, int minSupport);
std::vector<variant_t> clusterAndAnalyze(const std::string &bamFile, const std::vector<std::string> &chromosomes, const std::string &refFasta, const std::string &targetsBed, const std::string &excludeBed, int minOpSize, int maxSeparation, int minSupport, int numThreads, bool contigsOut);

int main(int argc, char *argv[]) {
    std::string tumorBamFile;
    std::string normalBamFile;
    std::string refFasta;
    std::string vcfNormal;
    std::string vcfTumor;
    std::string targetsBed;
    std::string excludeBed;

    int minOpSize = 1;
    int maxSeparation = 5;
    int minSupport = 3;
    int numThreads = 1;
    bool contigsOut = false;

    // Basic command-line argument parsing
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--tumor-bam" && i + 1 < argc) {
            tumorBamFile = argv[++i];
        } else if (arg == "--normal-bam" && i + 1 < argc) {
            normalBamFile = argv[++i];
        } else if (arg == "--ref" && i + 1 < argc) {
            refFasta = argv[++i];
        } else if (arg == "--bed" && i + 1 < argc) {
            targetsBed = argv[++i];
        } else if (arg == "--exclude-bed" && i + 1 < argc) {
            excludeBed = argv[++i];
        } else if (arg == "--normal-vcf" && i + 1 < argc) {
            vcfNormal = argv[++i];
        } else if (arg == "--tumor-vcf" && i + 1 < argc) {
            vcfTumor = argv[++i];
        } else if (arg == "--contigs-out" && i < argc) {
            contigsOut = true;
        } else if (arg == "--minOpSize" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> minOpSize;
        } else if (arg == "--maxSeparation" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> maxSeparation;
        } else if (arg == "--minSupport" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> minSupport;
        } else if (arg == "--threads" && i + 1 < argc) {
            std::istringstream(argv[++i]) >> numThreads;
        } else {
            std::cerr << " ERROR: Unknown argument or missing value for: " << arg << std::endl;
            return 1;
        }
    }

    if (tumorBamFile.empty() || refFasta.empty() || vcfTumor.empty()) {
        std::cerr << "Usage: " << "./rivindel" << std::endl << "Version:." << std::endl << std::endl
                  << "Input parameters" << std::endl
                  << "   --tumor-bam   STRING  Input tumor BAM file" << std::endl
                  << "   --normal-bam  STRING  Input normal BAM file (optional for paired analysis)" << std::endl
                  << "   --ref         STRING  Reference genome in FASTA format" << std::endl
                  << "   --tumor-vcf   STRING  Output tumor VCF file" << std::endl
                  << "   --normal-vcf  STRING  Output normal VCF file (optional for paired analysis)" << std::endl
                  << "   --threads     INTEGER Number of threads (default=1)" << std::endl << std::endl
                  << "For targeted sequencing:" << std::endl
                  << "   --bed         STRING  Targeted regions BED file" << std::endl << std::endl
                  << "Thresholds: " << std::endl
                  << "   --minOpSize   INTEGER Minimum number of CIGAR operations to consider a complex Indel" << std::endl
                  << "   --maxSeparation INTEGER Minimum separation in bp between CIGAR operations" << std::endl
                  << "   --minSupport  INTEGER Minimum read support" << std::endl;
        return 1;
    }

    // Parse target regions from BED file
    auto targetRegions = parseBedFile(targetsBed);

    std::set<std::string> targetChromosomes;
    for (const auto& region : targetRegions) {
        targetChromosomes.insert(region.chromosome);
    }

    // Get all chromosomes from the BAM file
    BamReader reader(tumorBamFile);
    auto header = reader.getHeader();
    std::vector<std::string> chromosomes;

    if (header) {
        for (int i = 0; i < header->n_targets; ++i) {
            std::string chrName = header->target_name[i];
            if (targetChromosomes.find(chrName) != targetChromosomes.end()) {
                chromosomes.push_back(chrName);
            }
        }
    }

    // Step 1: Extract signals for tumor BAM
    extractSignals(tumorBamFile, chromosomes, targetsBed, excludeBed, minOpSize, maxSeparation, minSupport);

    // Step 2: Cluster and analyze for tumor BAM
    std::vector<variant_t> tumorVariants = clusterAndAnalyze(tumorBamFile, chromosomes, refFasta, targetsBed, 
        excludeBed, minOpSize, maxSeparation, minSupport, numThreads, contigsOut);

    // Write all collected variants to the tumor VCF file
    PopulateVCF(tumorVariants, tumorBamFile, vcfTumor);

    // If a normal BAM file is provided, analyze it as well
    if (!normalBamFile.empty() && !vcfNormal.empty()) {
        // Extract signals for normal BAM
        extractSignals(normalBamFile, chromosomes, targetsBed, excludeBed, minOpSize, maxSeparation, minSupport);

        // Cluster and analyze for normal BAM
        std::vector<variant_t> normalVariants = clusterAndAnalyze(normalBamFile, chromosomes, refFasta, targetsBed, 
            excludeBed, minOpSize, maxSeparation, minSupport, numThreads, contigsOut);

        // Write all collected variants to the normal VCF file
        PopulateVCF(normalVariants, normalBamFile, vcfNormal);


        std::string outputVcf = vcfTumor + ".final.vcf";

       // Compare tumor and normal variants and classify them
        classififyVariants(vcfTumor, tumorBamFile, vcfNormal, outputVcf);
    }
    return 0;
}

void extractSignals(const std::string &bamFile, const std::vector<std::string> &chromosomes, const std::string &targetsBed, const std::string &excludeBed, int minOpSize, int maxSeparation, int minSupport) {
    int testThreads = 2;
    #pragma omp parallel for num_threads(testThreads) schedule(dynamic)
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        const auto& chromosomeName = chromosomes[i];

        std::string tempBamFile = bamFile + "." + chromosomeName + ".signals.bam";
        if (std::filesystem::exists(tempBamFile)) {
            continue;
        }

        // Create temporary BAM file for this chromosome
        createTempBamFile(bamFile, chromosomeName, tempBamFile, targetsBed, excludeBed, minOpSize, maxSeparation, minSupport);
    }
}

std::vector<variant_t> clusterAndAnalyze(const std::string &bamFile, const std::vector<std::string> &chromosomes, const std::string &refFasta, const std::string &targetsBed, const std::string &excludeBed, int minOpSize, int maxSeparation, int minSupport, int numThreads, bool contigsOut) {
    std::vector<variant_t> allVariants;
    std::mutex variantsMutex;

    // #pragma omp parallel for num_threads(numThreads) schedule(dynamic)
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        const auto& chromosomeName = chromosomes[i];
        std::string tempBamFile = bamFile + "." + chromosomeName + ".signals.bam";
        std::string clustersLog = bamFile + "." + chromosomeName + ".clusters.txt";

        // Cluster reads for this chromosome
        std::map<std::string, std::vector<clustered_aln_t>> clusteredReads = clusterInformativeReads(
            bamFile, tempBamFile, minOpSize, maxSeparation, minSupport, targetsBed, excludeBed, clustersLog, chromosomeName);

        std::cout << " INFO: Processing " << chromosomeName << std::endl;
        // Detect indels for this chromosome
        std::vector<variant_t> complexVariants = DetectComplexIndels(clusteredReads, bamFile, refFasta, numThreads, chromosomeName, contigsOut, "");

        // Append variants for this chromosome to the total list
        std::lock_guard<std::mutex> lock(variantsMutex);
        allVariants.insert(allVariants.end(), complexVariants.begin(), complexVariants.end());
    }

    return allVariants;
}

