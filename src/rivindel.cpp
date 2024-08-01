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

// #include "cluster_t.h"

int main(int argc, char *argv[]) {
    std::string bamFile;
    std::string refFasta;
    std::string vcfOut;
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
        if (arg == "--bam" && i + 1 < argc) {
            bamFile = argv[++i];
        } 
        else if (arg == "--ref" && i + 1 < argc) {
            refFasta = argv[++i];
        }
        else if (arg == "--bed" && i + 1 < argc) {
            targetsBed = argv[++i];
        }
        else if (arg == "--exclude-bed" && i + 1 < argc) {
            excludeBed = argv[++i];
        }
        else if (arg == "--vcf" && i + 1 < argc) {
            vcfOut = argv[++i];
        }
        else if (arg == "--contigs-out" && i < argc) {
            contigsOut = true;
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

    // Parse target regions from BED file
    auto targetRegions = parseBedFile(targetsBed);

    std::set<std::string> targetChromosomes;
    for (const auto& region : targetRegions) {
        targetChromosomes.insert(region.chromosome);
    }

    // Get all chromosomes from the BAM file
    BamReader reader(bamFile);
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
    std::string contigsOutFile;
    if (contigsOut) {

        contigsOutFile = bamFile + "." + "contigs.txt";
        std::ofstream outFile(contigsOutFile);
    }

    std::vector<variant_t> allVariants;
    for (const auto& chr : chromosomes) {
        std::string chromosomeName = chr;
        std::cout << " INFO: Clustering " << chromosomeName << std::endl;

        std::string tempBamFile = bamFile + "." + chromosomeName + ".signals.bam";
        std::string clustersLog = bamFile + "." + chromosomeName + ".clusters.txt";

        // Cluster reads for this chromosome
        std::map<std::string, std::vector<clustered_aln_t>> clusteredReads = clusterInformativeReads(
            bamFile, tempBamFile, minOpSize, maxSeparation, minSupport, targetsBed, excludeBed, clustersLog, chromosomeName);

        std::cout << " INFO: Processing " << chromosomeName << std::endl;

        // Detect indels for this chromosome
        std::vector<variant_t> complexVariants = DetectComplexIndels(clusteredReads, bamFile, refFasta, 
            numThreads, chromosomeName, contigsOut, contigsOutFile);

        // Append variants for this chromosome to the total list
        allVariants.insert(allVariants.end(), complexVariants.begin(), complexVariants.end());
    }

    // Write all collected variants to the VCF file
    PopulateVCF(allVariants, bamFile, vcfOut);

    return 0;
}


// int main(int argc, char *argv[]) {

//     std::string bamFile;
//     std::string refFasta;
//     std::string vcfOut;
//     std::string targetsBed;

//     int minOpSize = 1;
//     int maxSeparation = 5; 
//     int minSupport = 3;
//     int numThreads = 1;


//     // Basic command-line argument parsing
//     for (int i = 1; i < argc; ++i) {
//         std::string arg = argv[i];
//         if (arg == "--bam" && i + 1 < argc) {
//             bamFile = argv[++i];
//         } 
//         else if (arg == "--ref" && i + 1 < argc) {
//             refFasta = argv[++i];
//         }
//         else if (arg == "--bed" && i + 1 < argc) {
//             targetsBed = argv[++i];
//         }
//         else if (arg == "--vcf" && i + 1 < argc) {
//             vcfOut = argv[++i];
//         } 
//         else if (arg == "--minOpSize" && i + 1 < argc) {
//             std::istringstream(argv[++i]) >> minOpSize;
//         } 
//         else if (arg == "--maxSeparation" && i + 1 < argc) {
//             std::istringstream(argv[++i]) >> maxSeparation;
//         } 
//         else if (arg == "--minSupport" && i + 1 < argc) {
//             std::istringstream(argv[++i]) >> minSupport;
//         } 
//         else if (arg == "--threads" && i + 1 < argc) {
//             std::istringstream(argv[++i]) >> numThreads;
//         } 

//         else {
//             std::cerr << " ERROR: Unknown argument or missing value for: " << arg << std::endl;
//             return 1;
//         }
//     }

//     if (bamFile.empty() || refFasta.empty() || vcfOut.empty()) {
//         std::cerr << "Usage: " << "./rivindel" << std::endl << "Version:." << std::endl << std::endl <<
//                 "Input parameters" << std::endl <<
//                 "   --bam   STRING  Input BAM file" << std::endl <<
//                 "   --ref   STRING  Reference genome in FASTA format" << std::endl <<
//                 "   --vcf   STRING  Output VCF file" << std::endl <<
//                 "   --threads   INTEGER Number of threads (default=1)" << std::endl << std::endl <<

//                 "For targeted sequencing:" << std::endl <<
//                 "   --bed   STRING  Targeted regions BED file" << std::endl << std::endl <<

//                 "Thresholds: " << std::endl <<
//                 "   --minOpSize INTEGER Minimum number of CIGAR operations to consider a complex Indel" << std::endl <<
//                 "   --maxSeparation INTEGER Minimum separation in bp between CIGAR operations" << std::endl <<
//                 "   --minSupport    INTEGER Minimum read support" << std::endl;
//         return 1;
//     }

//     std::filesystem::path bamPath(bamFile);
//     std::string baseName = bamPath.stem().string(); // Get the basename without extension
//     std::string tempBamFile = bamPath.parent_path().string() + "/" + baseName + ".rivindel.signals.bam"; // Appe

//     std::string clustersLog = bamPath.parent_path().string() + "/" + baseName + ".rivindel.clusters.txt"; // Appe


//     // identify cluster reads with compatible signals for cxindels
//     std::map<std::string, std::vector<clustered_aln_t>> clusteredReads = clusterInformativeReads(bamFile, 
//         minOpSize, maxSeparation, minSupport, tempBamFile, targetsBed, clustersLog);


//     // Assemble clusters, align to reference and call variants
//     std::vector<variant_t> complexVariants = DetectComplexIndels(clusteredReads, 
//         bamFile, refFasta, numThreads);

//     PopulateVCF(complexVariants, bamFile, vcfOut);

//     return 0;
// }