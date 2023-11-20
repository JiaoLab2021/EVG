#include "../include/check_vcf_sort.hpp"

using namespace std;

int check_vcf_sort(string inputVcf)
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Checking sort: '"<< inputVcf << "'\n";

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(inputVcf);

    // Record the start and chromosome number of the last mutation
    string preChromosome;
    uint32_t preRefStart = 0;

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty or comment line, skip
        if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        // If it is a new chromosome, zero the starting position
        if (INFOSTRUCTTMP.CHROM != preChromosome) {
            preChromosome = INFOSTRUCTTMP.CHROM;
            preRefStart = 0;
        }

        // greater than or equal to the previous variation
        if (INFOSTRUCTTMP.POS >= preRefStart) {
            preRefStart = INFOSTRUCTTMP.POS;
        } else {  // Smaller than the previous variation, it means that the vcf is not sorted, and the user is reminded to sort
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: unsorted file -> "<< inputVcf << " " 
                << INFOSTRUCTTMP.CHROM << " " << preRefStart << ">" << INFOSTRUCTTMP.POS << endl;
            exit(1);
        }
    }

    return 0;
}