// g++ vcf_convert.cpp -o vcf_convert -lz

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv) {
    // Enter the file and parameters
    string refFileName;
    string vcfFileName;
    string outFileName = "vcfConvert.out.vcf.gz";
    int readLen = 350; // read length

    // Filtering threshold
    double MAF = 0;  // Minimum allele frequency
    double MISSRATE = 1.0;  // Miss rate

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},

            {"length", required_argument, 0, 'l'},
            {"maf", required_argument, 0, '1'},
            {"geno", required_argument, 0, '2'},
            {"out", required_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "r:v:l:1:2:o:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            refFileName = optarg;
            break;
        case 'v':
            vcfFileName = optarg;
            break;
        case 'l':
            readLen = stoi(optarg);
            break;
        case '1':
            MAF = stod(optarg);
            break;
        case '2':
            MISSRATE = stod(optarg);
            break;
        case 'o':
            outFileName = optarg;
            break;
        case 'h':
        case '?':
            help_convert(argv);
            exit(1);
            break;
        default:
            abort();
        }
    }

    if (argc <= 2) {
        help_convert(argv);
        return 1;
    }

    // Print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // If the length of "raed" is greater than 1000, there's no need to filter by read length because we don't have to process paragraph
    readLen = (readLen > 1000) ? 0 : readLen;

    // init
    Convert ConvertClass(refFileName, vcfFileName, readLen, MAF, MISSRATE, outFileName);
    // build reference index
    ConvertClass.build_reference_index();
    // Convert VCF file
    ConvertClass.vcf_convert();

    // Print log
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// Help document
void help_convert(char** argv) {
  cerr << "usage: " << argv[0] << " " << argv[1] << " -r FILE -v FILE [options]" << endl
       << "convert VCF files to the required format for genome graph software." << endl
       << endl
       << "required arguments:" << endl
       << "    -r, --reference   FILE     input FASTA reference" << endl
       << "    -v, --vcf         FILE     input VCF" << endl
       << endl
	   << "optional arguments:" << endl
       << "    -l, --length      INT      read length [350]" << endl
       << "    --maf             FLOAT    exclude SNPs with minor allele frequency lower than threshold [0.0]" << endl
       << "    --geno            FLOAT    exclude SNPs with missing call frequencies greater than threshold [1.0]" << endl
       << "    -o, --out         FILE     output file name [vcfConvert.out.vcf.gz]" << endl
       << endl
       << "    -h, --help                 print this help document" << endl;
}



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
Convert::Convert(
    string refFileName, 
    string vcfFileName, 
    int readLen, 
    double MAF, 
    double MISSRATE, 
    string outFileName
) : refFileName_(refFileName), vcfFileName_(vcfFileName), readLen_(readLen), MAF_(MAF), MISSRATE_(MISSRATE), outFileName_(outFileName) {}


/**
 * @brief build the reference genome index
 * 
 * @return void
**/
void Convert::build_reference_index() {
    // Chromosome name output, graphtyper required
    ofstream outFile;
    outFile.open("CHROMOSOME.NAME", ios::out);
    if(!outFile.is_open()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'CHROMOSOME.NAME': No such file or directory." 
            << endl;
        exit(1);
    }

    // Input file stream
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // open file
    if(!gzfp) {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    } else {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while(kseq_read(ks) >= 0) {
            string chromosome = ks->name.s;
            uint32_t chrLen = ks->seq.l;
            string sequence = ks->seq.s;

            refIndexS_.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
            refIndexS_.chrSeqMap[chromosome] = sequence;
            refIndexS_.chrLenMap[chromosome] = chrLen;

            // Output chromosome name
            outFile << chromosome + "\n";
        }

        // Free the memory and close the file
        kseq_destroy(ks);
        gzclose(gzfp);
    }

    // free memory
    outFile.close();
}


/**
 * 0. Check Line number consistency
 * 1. Head plus chromosome length. (bayestyper)
 * 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
 * 3. The first variation is greater than the read length. (paragraph)
 * 4. The sequence must correspond to reference.
 * 5. Replace 'END=' in the eighth column of vcf. (paragrpah)
 * 6. Graphtyper2 will report an error when it encounters SVLEN=1,2.
 * 7. The first base of refSeq is the same as the first base of qrySeq.
 * 8. Check whether ref and qry sequences are the same and the same sites are skipped.
 * 9. Check whether refSeq and qrySeq contain characters other than atgcnATGCN. If yes, skip this site.
 * 10. Check whether the qry after replacement is the same. If yes, skip this site. (BayesTyper, Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed)
 * 11. Combine the genotypes. Convert to.|. (PanGenie)
 * 12. Change the / in genotype to |. (PanGenie)
 * 13. Retain only the diploid variation in the genotype.
 * 14. Check if GT has more sequences than qry.
 * 15. If the variant end position is greater than or equal to the chromosome length, skip the mutation. (Paragrpah)
 * 16. Duplicate positions.
 * 17. Trim ALT alleles not seen in the genotype fields and update the genotypes accordingly.
**/
/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @return void
**/
void Convert::vcf_convert() {
    uint32_t gtIndex = 0;  // Genotype index

    // Regular expression
    const std::regex endReg(R"(END=\d+)");
    const std::regex svlenReg1(R"(SVLEN=[-\d+,]*;)");
    const std::regex svlenReg2(R"(;SVLEN=[-\d+,]*)");
    const std::regex svlenReg3(R"(SVLEN=[-\d+,]*)");

    // Sequence information on the genome
    const map<string, string>& chrSeqMap = refIndexS_.chrSeqMap;
    const map<string, uint32_t>& chrLenMap = refIndexS_.chrLenMap;

    // Output file stream
    SAVE SAVECLASS(outFileName_);

    stringstream outStream; // Use stringstream instead of string concatenation
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // Input file stream
    // Store vcf information
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    uint32_t lineNum = 0;  // Line number

    // Record the information of the previous variant
    VCFINFOSTRUCT prePreINFOSTRUCTTMP;
    VCFINFOSTRUCT preINFOSTRUCTTMP;

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        uint32_t thisLineNum = INFOSTRUCTTMP.lineVec.size();
        
        // comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
                // 1. Head plus chromosome length. (bayestyper)
                outStream << refIndexS_.chrLenTxt;

                // save header
                outStream << INFOSTRUCTTMP.line + "\n";

                // Line number
                lineNum = thisLineNum;
            }
            // Skip if the vcf contains chromosome length information
            else if (INFOSTRUCTTMP.line.find(",length") == string::npos) {
                outStream << INFOSTRUCTTMP.line + "\n";
            }

            continue;
        }


        // 0. Check Line number consistency
        if (lineNum > 0 && lineNum != thisLineNum) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: The number of columns in the VCF file is inconsistent. Skipping this site -> " 
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // Find the location of gt in the FORMAT field
        vector<string> formatVec = split(INFOSTRUCTTMP.FORMAT, ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        // GT index
        gtIndex = distance(formatVec.begin(), gtItera);
        if (gtIndex == formatVec.size()) {  // Check whether index exists. If no index exists, the genotypes are all 0.
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "Warning: The FORMAT field does not contain genotype (GT) information at the following location: " 
                 << INFOSTRUCTTMP.CHROM << ":" << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        /* ************************ filter SNPs by maf and missing rate ************************ */
        // Determine whether to filter
        if (MAF_ > 0 && MISSRATE_ < 1) {
            // Acquire variant type
            INFOSTRUCTTMP.ID = VCFOPENCLASS.get_TYPE(
                INFOSTRUCTTMP.LEN, 
                INFOSTRUCTTMP.ALTVec
            );

            if (INFOSTRUCTTMP.ID == "SNP") {  // Determine whether to filter when SNP is displayed
                double MAFTMP;  // Minimum allele frequency
                double MISSRATETMP;  // Miss rate

                // Get all genotypes map<idx, vector<gtString> >
                unordered_map<int, vector<string> > GTVecMapTmp;
                uint32_t misSampleNum;
                tie(GTVecMapTmp, misSampleNum) = VCFOPENCLASS.get_gt(INFOSTRUCTTMP.lineVec);

                // If there is only one genotype, skip.
                if (GTVecMapTmp.size() - misSampleNum <= 1) continue;
                
                // The minimum allele frequency and deletion rate were calculated
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    misSampleNum
                );

                // Did not pass the threshold, straight to the next loop
                if (MAFTMP < MAF_ || MISSRATETMP > MISSRATE_) continue;
            }
        }


        // Check if there is corresponding chromosome information in the submitted genome
        if (chrSeqMap.find(INFOSTRUCTTMP.CHROM) == chrSeqMap.end() || chrLenMap.find(INFOSTRUCTTMP.CHROM) == chrLenMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The chromosome '"
                << INFOSTRUCTTMP.CHROM 
                << "' is not found in the reference file '" << refFileName_ << "'" 
                << endl;
            continue;
        }

        // 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
        if (INFOSTRUCTTMP.lineVec[4].find("<") != string::npos) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The query sequence contains the symbol '<'. Skipping this site -> " 
                << INFOSTRUCTTMP.CHROM << " " 
                << INFOSTRUCTTMP.POS << " " 
                << INFOSTRUCTTMP.lineVec[4] << endl;
            continue;
        }
        

        // 3. The first variation is greater than the read length. (paragraph)
        if (static_cast<uint32_t>(readLen_) > INFOSTRUCTTMP.POS) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The variant's start position is less than the read length. Skipping this site -> "
                << INFOSTRUCTTMP.CHROM << " " 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // 4. The sequence must correspond to reference.
        if (chrLenMap.at(INFOSTRUCTTMP.CHROM) < (INFOSTRUCTTMP.END)) {  // Check that the chromosome length is correct
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The variant's end position exceeds the length of the chromosome -> " 
                << INFOSTRUCTTMP.CHROM << " " 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        string trueRefSeq = chrSeqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) {  // If it is different from the sequence in the reference genome, it is replaced with the sequence in the reference genome
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: Sequence discrepancy between the reference genome and the VCF, replacing with the sequence from the reference genome -> "
                << INFOSTRUCTTMP.CHROM << " " 
                << INFOSTRUCTTMP.POS << endl;

            INFOSTRUCTTMP.lineVec[3] = trueRefSeq;
            INFOSTRUCTTMP.REF = trueRefSeq;
        }


        // 5. Replace the 'END=' in the eighth column of the VCF. (paragrpah)
        // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end))) ("{}:{} missing REF or ALT sequence.".format(start, end))
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 6. Graphtyper2 will report an error when it encounters SVLEN=1,2.
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg1, "");
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg2, "");
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg3, ".");


        // 7. The first base of refSeq is the same as the first base of qrySeq. (paragraph)
        // 8. Check whether ref and qry sequences are the same and the same sites are skipped
        // 9. Check whether refSeq and qrySeq contain characters other than atgcnATGCN. If yes, skip this site.
        process_sequences(INFOSTRUCTTMP, chrSeqMap);


        // 10. Check whether the qry after replacement is the same. If yes, skip this site. (BayesTyper, Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed)
        check_repeated_alleles(INFOSTRUCTTMP);


        uint32_t maxGT = 0; // Record the largest GT, see if there are more sequences than qry, and skip this site if there are more
        // As long as GT field
        INFOSTRUCTTMP.lineVec[8] = "GT";
        if (gtItera != formatVec.end()) {  // The FORMAT contains GT
            // Loop over the GT column
            for (size_t i = 9; i < thisLineNum; i++) {
                string gt;

                string genotypeInfo = INFOSTRUCTTMP.lineVec[i];
                if (genotypeInfo.empty()) {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: Genotype (GT) is empty ->" 
                        << INFOSTRUCTTMP.CHROM << " " 
                        << INFOSTRUCTTMP.POS << endl;
                    gt = "0|0";
                } else {
                    // Find the genotype
                    gt = split(genotypeInfo, ":")[gtIndex];
                }

                // 11. Combine the genotypes. Convert to.|.. (PanGenie)
                if (gt == ".") {
                    gt = ".|.";
                } else if (gt == "0") {
                    gt = "0|0";
                }


                // 12. Change the / in genotype to |. (PanGenie)
                if (gt.find("/") != string::npos) {
                    std::regex reg("/");
                    gt = regex_replace(string(gt), regex(reg), string("|"));
                }

                // Loop for the largest GT
                vector<string> gtVec = VCFOPENCLASS.gt_split(gt);
                for (auto it1 : gtVec) {
                    try {
                        if (stoul(it1) > maxGT) {
                            maxGT = stoi(it1);
                        }
                    } catch(const std::invalid_argument& e) {
                        continue;
                    }
                }

                // 13. Retain only the diploid variation in the genotype
                if (gtVec.size() == 1) {
                    gt = gtVec[0] + "|0";
                } else if (gtVec.size() > 2) {
                    gt = gtVec[0] + "|" + gtVec[1];
                }
                
                INFOSTRUCTTMP.lineVec[i] = gt;
            }
        } else {  // Skip the site
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: Genotype (GT) not found in FORMAT column. Skipping this site -> " 
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 14. Check if GT has more sequences than qry
        // pangenie: VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: The number of genotypes (GTs) exceeds the number of query sequences at: " 
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 15. If the variant end position is greater than or equal to the chromosome length, skip the mutation. (Paragrpah)
        if (chrLenMap.at(INFOSTRUCTTMP.CHROM) == (INFOSTRUCTTMP.END)) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: The variant's end position equals the chromosome length. Skipping this site -> " 
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // Check whether the information is cleared. If it is cleared, skip it
        if (INFOSTRUCTTMP.line.empty()) continue;


        // 16. Duplicate positions
        if (INFOSTRUCTTMP.CHROM == preINFOSTRUCTTMP.CHROM && INFOSTRUCTTMP.POS == preINFOSTRUCTTMP.POS) {
            // Update genotype information
            merge_loci(
                gtIndex, 
                preINFOSTRUCTTMP, 
                INFOSTRUCTTMP, 
                VCFOPENCLASS, 
                chrSeqMap
            );
            continue;
        } else if (INFOSTRUCTTMP.CHROM == prePreINFOSTRUCTTMP.CHROM && INFOSTRUCTTMP.POS == prePreINFOSTRUCTTMP.POS) {
            // Update genotype information
            merge_loci(
                gtIndex,
                prePreINFOSTRUCTTMP,
                INFOSTRUCTTMP,
                VCFOPENCLASS, 
                chrSeqMap
            );
            continue;
        } else if (INFOSTRUCTTMP.CHROM == prePreINFOSTRUCTTMP.CHROM && INFOSTRUCTTMP.POS > prePreINFOSTRUCTTMP.POS && INFOSTRUCTTMP.POS < preINFOSTRUCTTMP.POS) {
            // If the current node is smaller than the previous node after replacing the coordinates, but larger than the previous node, replace the current node with the previous node.
            INFOSTRUCTTMP.swap(preINFOSTRUCTTMP);
        } else if (INFOSTRUCTTMP.CHROM == prePreINFOSTRUCTTMP.CHROM && INFOSTRUCTTMP.POS < prePreINFOSTRUCTTMP.POS) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The variant's start position is less than the previous one. Skipping this site -> " 
                << INFOSTRUCTTMP.CHROM << " " 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // Add the replaced string to outStream
        if (INFOSTRUCTTMP.CHROM != prePreINFOSTRUCTTMP.CHROM || prePreINFOSTRUCTTMP.POS != INFOSTRUCTTMP.POS) {
            if (!prePreINFOSTRUCTTMP.CHROM.empty()) {
                // 17. Trim ALT alleles not seen in the genotype fields and update the genotypes accordingly.
                trim_update_alt_alleles(VCFOPENCLASS, prePreINFOSTRUCTTMP, gtIndex);

                // 10. Check whether the qry after replacement is the same. If yes, skip this site. (BayesTyper, Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed)
                check_repeated_alleles(prePreINFOSTRUCTTMP);

                if (!prePreINFOSTRUCTTMP.line.empty() && !prePreINFOSTRUCTTMP.ALTVec.empty()) {  // (BayesTyper, Assertion `!alt_allele->empty()' failed)
                    outStream << join(prePreINFOSTRUCTTMP.lineVec, "\t") + "\n";
                }
            }
            prePreINFOSTRUCTTMP = std::move(preINFOSTRUCTTMP);
            preINFOSTRUCTTMP = std::move(INFOSTRUCTTMP);
        }

        if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
            string outTxt = outStream.str();
            SAVECLASS.save(outTxt);
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    // 17. Trim ALT alleles not seen in the genotype fields and update the genotypes accordingly.
    trim_update_alt_alleles(VCFOPENCLASS, prePreINFOSTRUCTTMP, gtIndex);
    trim_update_alt_alleles(VCFOPENCLASS, preINFOSTRUCTTMP, gtIndex);

    // 10. Check whether the qry after replacement is the same. If yes, skip this site. (BayesTyper, Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed)
    check_repeated_alleles(prePreINFOSTRUCTTMP);
    check_repeated_alleles(preINFOSTRUCTTMP);

    if (!prePreINFOSTRUCTTMP.line.empty() && !prePreINFOSTRUCTTMP.ALTVec.empty()) {
        outStream << join(prePreINFOSTRUCTTMP.lineVec, "\t") + "\n";
    }
    if (!preINFOSTRUCTTMP.line.empty() && !preINFOSTRUCTTMP.ALTVec.empty()) {
        outStream << join(preINFOSTRUCTTMP.lineVec, "\t") + "\n";
    }

    if (outStream.tellp() > 0) {  //Write for the last time
        string outTxt = outStream.str();
        SAVECLASS.save(outTxt);
        // Clear stringstream
        outStream.str(string());
        outStream.clear();
    }
}


/**
 * Update genotype information
 * 
 * @param gtIndex             Genotype index
 * @param preINFOSTRUCTTMP    Previous site information
 * @param INFOSTRUCTTMP       Current location information
 * @param VCFOPENCLASS        
 * @param chrSeqMap           A map from chromosome names to sequences.
 * 
 * @return void
*/
void Convert::merge_loci(
    const uint32_t& gtIndex, 
    VCFINFOSTRUCT& preINFOSTRUCTTMP, 
    VCFINFOSTRUCT& INFOSTRUCTTMP, 
    VCFOPEN& VCFOPENCLASS, 
    const std::map<std::string, std::string>& chrSeqMap
) {
    // if (std::abs(static_cast<int>(preINFOSTRUCTTMP.REF.size()) - static_cast<int>(INFOSTRUCTTMP.REF.size())) > 50) {
    //     cerr << "[" << __func__ << "::" << getTime() << "] "
    //     << "Warning: multiple variants observed, skip this site -> " 
    //     << INFOSTRUCTTMP.CHROM << " " 
    //     << INFOSTRUCTTMP.POS << endl;
    //     return;
    // }

    int32_t preAltNum = preINFOSTRUCTTMP.ALTVec.size();

    bool isPreLonger = preINFOSTRUCTTMP.REF.size() > INFOSTRUCTTMP.REF.size();
    string REF = isPreLonger ? preINFOSTRUCTTMP.REF : INFOSTRUCTTMP.REF;
    string diffREF = REF.substr(std::min(preINFOSTRUCTTMP.REF.size(), INFOSTRUCTTMP.REF.size()));

    // Get all genotypes map<idx, vector<gtString> >
    unordered_map<int, vector<string> > preGTVecMapTmp, GTVecMapTmp;
    uint32_t misSampleNum;
    tie(preGTVecMapTmp, misSampleNum) = VCFOPENCLASS.get_gt(preINFOSTRUCTTMP.lineVec);
    tie(GTVecMapTmp, misSampleNum) = VCFOPENCLASS.get_gt(INFOSTRUCTTMP.lineVec);

    // replace genotype
    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        vector<string> preGtVec = preGTVecMapTmp.at(i - 9);
        vector<string> gtVec = GTVecMapTmp.at(i - 9);

        // sort + uniqe + remove
        vector<string> preGtVecTmp = sort_unique(preGtVec);
        vector<string> gtVecTmp = sort_unique(gtVec);

        // If the previous genotype is not 0 or this genotype is 0, skip
        if ((preGtVecTmp.back() != "0" && preGtVecTmp.back() != ".") || (gtVecTmp.back() == "0" || gtVecTmp.back() == ".")) {
            continue;
        }

        // replace GT
        for (size_t j = 0; j < preGtVec.size(); j++) {
            if (gtVec[j] == "0" || gtVec[j] == "."){
                continue;
            }
            preGtVec[j] = to_string(stoul(gtVec[j]) + preAltNum);
        }
        vector<string> sampleVec = split(preINFOSTRUCTTMP.lineVec[i], ":");
        sampleVec[gtIndex] = join(preGtVec, "|");
        preINFOSTRUCTTMP.lineVec[i] = join(sampleVec, ":");
    }

    // replace REF
    preINFOSTRUCTTMP.REF = REF;
    preINFOSTRUCTTMP.lineVec[3] = REF;

    // replace LEN
    preINFOSTRUCTTMP.LEN = preINFOSTRUCTTMP.REF.size();

    // replace ALT
    if (isPreLonger) {
        for (auto & altSeq : INFOSTRUCTTMP.ALTVec) {
            altSeq += diffREF;
        }
    } else {
        for (auto & altSeq : preINFOSTRUCTTMP.ALTVec) {
            altSeq += diffREF;
        }
    }
    preINFOSTRUCTTMP.ALTVec.insert(preINFOSTRUCTTMP.ALTVec.end(), INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end());
    preINFOSTRUCTTMP.ALT = join(preINFOSTRUCTTMP.ALTVec, ",");
    preINFOSTRUCTTMP.lineVec[4] = preINFOSTRUCTTMP.ALT;

    // replace END
    preINFOSTRUCTTMP.END = preINFOSTRUCTTMP.POS + preINFOSTRUCTTMP.REF.size() - 1;

    // replace INFO
    static const std::regex endReg(R"(END=\d+)");  // Regular expression
    preINFOSTRUCTTMP.lineVec[7] = regex_replace(preINFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(preINFOSTRUCTTMP.END));

    // 7. The first base of refSeq is the same as the first base of qrySeq. (paragraph)
    process_sequences(preINFOSTRUCTTMP, chrSeqMap);
}

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
void Convert::process_sequences(VCFINFOSTRUCT& INFOSTRUCTTMP, const std::map<std::string, std::string>& chrSeqMap) {
    for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
        string qrySeq = INFOSTRUCTTMP.ALTVec[i];

        // 7. The first base of refSeq is the same as the first base of qrySeq. (paragraph)
        if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1)) {
            // Move the coordinates one step backward
            INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
            INFOSTRUCTTMP.lineVec[1] = to_string(INFOSTRUCTTMP.POS);
            INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
            INFOSTRUCTTMP.REF = chrSeqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
            INFOSTRUCTTMP.lineVec[3] = INFOSTRUCTTMP.REF;  // Assigns a value to the Vector

            // Add 'refSeq[0]' to all sequences in 'INFOSTRUCTTMP.ALTVec'
            for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++) {
                INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
            }

            qrySeq = INFOSTRUCTTMP.ALTVec[i];
        }

        // 8. Check whether ref and qry sequences are the same and the same sites are skipped
        if (qrySeq == INFOSTRUCTTMP.REF) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: The REF and ALT sequences are the same. Skipping this site -> "
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            INFOSTRUCTTMP.line.clear();
            break;
        }

        // 9. Check whether refSeq and qrySeq contain characters other than atgcnATGCN. If yes, skip this site.
        smatch results;
        std::regex atgcReg("[^ATGCNatgcn]");
        if (regex_search(qrySeq, results, atgcReg)) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "Warning: Sequence contains non-ATGCNatgcn characters, skip this site -> "
                 << INFOSTRUCTTMP.CHROM << " " 
                 << INFOSTRUCTTMP.POS << endl;
            INFOSTRUCTTMP.line.clear();
            break;
        }
    }

    INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");
}


/**
 * @brief Trim ALT alleles not seen in the genotype fields and update the genotypes accordingly.
 * 
 * @param VCFOPENCLASS The VCFOPEN object containing the VCF data.
 * @param INFOSTRUCTTMP The INFOSTRUCTTMP object containing the VCF data.
 * @param gtIndex The index of the GT field in the FORMAT column.
 */
void Convert::trim_update_alt_alleles(VCFOPEN& VCFOPENCLASS, VCFINFOSTRUCT& INFOSTRUCTTMP, const uint32_t& gtIndex) {
    // Create a set to store the alleles that appear in the genotypes
    set<uint32_t> allelesInGenotypes;
    allelesInGenotypes.insert(0);  // Add the reference allele

    // Iterate over the samples
    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        // Split the genotype field
        string gt = split(INFOSTRUCTTMP.lineVec[i], ":")[gtIndex];

        // Split the alleles
        vector<string> alleles = VCFOPENCLASS.gt_split(gt);

        // Add the alleles to the set
        for (const string& allele : alleles) {
            try {
                allelesInGenotypes.insert(stoul(allele));
            } catch (const std::invalid_argument& e) {
                // Ignore non-integer alleles
            }
        }
    }

    // If the number of alleles in the genotypes is equal to the number of alleles (REF + ALT), there is no need to trim
    if (allelesInGenotypes.size() == INFOSTRUCTTMP.ALTVec.size() + 1) {
        return;
    }

    // Remove the alleles that do not appear in any genotype and update the genotypes
    vector<string> newALTVec;
    set<uint32_t> allelesToRemove;
    for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
        uint32_t alleleIndex = i + 1;  // +1 because REF is 0
        if (allelesInGenotypes.count(alleleIndex) == 0) {
            allelesToRemove.insert(alleleIndex);
        } else {
            newALTVec.push_back(INFOSTRUCTTMP.ALTVec[i]);
        }
    }

    // Update the ALT field
    INFOSTRUCTTMP.ALTVec = newALTVec;
    // Update the ALT field in the line vector
    INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");

    // Update the genotypes
    for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
        // Split the genotype field
        vector<string> gtVec = split(INFOSTRUCTTMP.lineVec[i], ":");
        string gt = gtVec[gtIndex];

        // Split the alleles
        vector<string> alleles = VCFOPENCLASS.gt_split(gt);

        // Update the alleles
        for (string& allele : alleles) {
            try {
                uint32_t oldAlleleNumber = stoul(allele);
                uint32_t newAlleleNumber = oldAlleleNumber;

                // For each allele to remove, if it is less than the current allele, decrement the new allele number
                for (const auto& alleleToRemove : allelesToRemove) {
                    if (alleleToRemove < oldAlleleNumber) {
                        --newAlleleNumber;
                    }
                }

                // Update the allele
                allele = to_string(newAlleleNumber);
            } catch (const std::invalid_argument& e) {
                // Ignore non-integer alleles
            }
        }

        // Update the genotype field
        gtVec[gtIndex] = join(alleles, "|");
        INFOSTRUCTTMP.lineVec[i] = join(gtVec, ":");
    }
}


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
void Convert::check_repeated_alleles(VCFINFOSTRUCT& infoStruct) {
    for (const auto& it : infoStruct.ALTVec) {
        if (count(infoStruct.ALTVec.begin(), infoStruct.ALTVec.end(), it) > 1) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: Repeated alleles found. Skipping this site -> "
                << infoStruct.CHROM << " " 
                << infoStruct.POS << endl;
            infoStruct.line.clear();
            break;
        }
    }
}


/**
 * @brief Sort, deduplicate and delete the vector
 * 
 * @param vec    
 * 
 * @return vec
*/
vector<string> Convert::sort_unique(vector<string> vec) {
    sort(vec.begin(), vec.end(), [](const string& a, const string& b) {
        return stoul(a) < stoul(b);
    });
    auto last = unique(vec.begin(), vec.end());
    vec.erase(last, vec.end());
    return vec;
}
