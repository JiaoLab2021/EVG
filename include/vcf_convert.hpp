#ifndef VCF_CONVERT_HPP
#define VCF_CONVERT_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include "zlib.h"
#include <getopt.h>
#include <sstream>
#include <unordered_set>
#include <set>

#include "strip_split_join.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "save.hpp"
#include "vcf_open.hpp"

using namespace std;

void help_convert(char** argv);

int main_convert(int argc, char** argv);


struct refIndexStruct {
    map<string, string> chrSeqMap;  // map<chromosome, sequence>
    map<string, uint32_t> chrLenMap;  // map<chromosome, length>
    string chrLenTxt;  // Store the contig length and save it to the vcf header

    refIndexStruct() : chrLenTxt("") {}
};


class Convert
{
private:
    string refFileName_;  // reference genome
    string vcfFileName_;  // input VCF file name

    int readLen_;  // read length

    double MAF_;  // Minimal allele frequency
    double MISSRATE_;  // Missing rate

    refIndexStruct refIndexS_;  // the index of reference genome

    string outFileName_;  // output file name

    // open file
    KSEQ_INIT(gzFile, gzread)

public:
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @param refFileName     reference genome
     * @param vcfFileName     input VCF file name
     * @param readLen         read length
     * @param MAF             Minimal allele frequency
     * @param MISSRATE        Missing rate
     * @param outFileName     output file name
     * 
    **/
    Convert(string refFileName, string vcfFileName, int readLen, double MAF, double MISSRATE, string outFileName);

    /**
     * @brief build the reference genome index
     * 
     * @return void
    **/
    void build_reference_index();

    /**
     * 1. Head plus chromosome length. (bayestyper)
     * 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
     * 3. The first variation is greater than the read length. (paragraph)
     * 4. The sequence must correspond to reference.
     * 5. Replace 'END=' in the eighth column of vcf. (paragrpah)
     * 6. Graphtyper2 will report an error when it encounters SVLEN=1,2.
     * 7. The first base of refSeq is the same as the first base of qrySeq.
     * 8. Check whether ref and qry sequences are the same and the same sites are skipped.
     * 9. Check whether refSeq and qrySeq contain characters other than atgcnATGCN. If yes, skip this site.
     * 10. Check whether the qry after replacement is the same. If yes, skip this site.
     * 11. Check for mutations that duplicate positions.
     * 12. Combine the genotypes. Convert to.|.. (PanGenie)
     * 13. Change the/in genotype to |. (PanGenie)
     * 14. Retain only the diploid variation in the genotype.
     * 15. Check if GT has more sequences than qry.
     * 16. If the variant end position is greater than or equal to the chromosome length, skip the mutation. (Paragrpah)
    **/
    /**
     * @brief Convert vcf to the format required by graph genome tools
     * 
     * @return void
    **/
    void vcf_convert();

    /**
     * Update genotype information
     * 
     * @param gtIndex             Genotype index
     * @param preINFOSTRUCTTMP    Previous site information
     * @param INFOSTRUCTTMP       Current location information
     * @param VCFOPENCLASS        
     * @param chrSeqMap A map from chromosome names to sequences.
     * 
     * @return void
    */
    void merge_loci(
        const uint32_t& gtIndex, 
        VCFINFOSTRUCT& preINFOSTRUCTTMP, 
        VCFINFOSTRUCT& INFOSTRUCTTMP, 
        VCFOPEN& VCFOPENCLASS, 
        const std::map<std::string, std::string>& chrSeqMap
    );

    /**
     * @brief Process the sequences in the given VCFINFOSTRUCT.
     * 
     * This function processes the sequences in the given VCFINFOSTRUCT. It adjusts the position and length of the sequences,
     * replaces the reference sequence if necessary, and checks for identical reference and alternate sequences. It also
     * checks for non-ATGCNatgcn characters in the sequences. If any of these checks fail, it clears the line and breaks
     * the loop.
     * 
     * @param INFOSTRUCTTMP The VCFINFOSTRUCT to process.
     * @param chrSeqMap A map from chromosome names to sequences.
    */
    void process_sequences(VCFINFOSTRUCT& INFOSTRUCTTMP, const std::map<std::string, std::string>& chrSeqMap);

    /**
     * @brief Trim ALT alleles not seen in the genotype fields and update the genotypes accordingly.
     * 
     * @param VCFOPENCLASS The VCFOPEN object containing the VCF data.
     * @param INFOSTRUCTTMP The INFOSTRUCTTMP object containing the VCF data.
     * @param gtIndex The index of the GT field in the FORMAT column.
    */
    void trim_update_alt_alleles(VCFOPEN& VCFOPENCLASS, VCFINFOSTRUCT& INFOSTRUCTTMP, const uint32_t& gtIndex);

    /**
     * @brief Checks for repeated alleles in the given INFOSTRUCT object and clears the line if any are found.
     * 
     * @param infoStruct An INFOSTRUCT object that contains the alleles to be checked.
     * 
     * This function iterates over the ALTVec member of the given INFOSTRUCT object. If it finds any repeated alleles,
     * it prints a warning message to the standard error stream, clears the line member of the INFOSTRUCT object, and
     * then breaks out of the loop.
     * 
     * Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed (BayesTyper)
     */
    void check_repeated_alleles(VCFINFOSTRUCT& infoStruct);

    /**
     * @brief Sort, deduplicate and delete the vector
     * 
     * @param vec    
     * 
     * @return vec
    */
    vector<string> sort_unique(vector<string> vec);
};

#endif