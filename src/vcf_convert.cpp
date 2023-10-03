// g++ vcf_convert.cpp -o vcf_convert -lz

#include "../include/vcf_convert.hpp"

using namespace std;


int main_convert(int argc, char** argv)
{
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
void help_convert(char** argv)
{
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
void Convert::build_reference_index()
{
    // Chromosome name output, graphtyper required
    ofstream outFile;
    outFile.open("CHROMOSOME.NAME", ios::out);
    if(!outFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'CHROMOSOME.NAME': No such file or directory." 
            << endl;
        exit(1);
    }

    // Input file stream
    gzFile gzfp = gzopen(refFileName_.c_str(), "rb");

    // open file
    if(!gzfp)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] "
            << "'" << refFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
    else
    {
        kseq_t *ks;
        ks = kseq_init(gzfp);
    
        while( kseq_read(ks) >= 0 )
        {
            string chromosome = ks->name.s;
            uint32_t chrLen = ks->seq.l;
            string sequence = ks->seq.s;

            refIndexS_.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
            refIndexS_.sequenceMap[chromosome] = sequence;

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
**/
/**
 * @brief Convert vcf to the format required by graph genome tools
 * 
 * @return void
**/
void Convert::vcf_convert()
{
    // Regular expression
    const std::regex endReg(R"(END=\d+)");
    const std::regex svlenReg(R"(SVLEN=[-\d+,]*;)");

    // Whether check file is sorted
    check_vcf_sort(vcfFileName_);

    // Sequence information on the genome
    const map<string, string>& seqMap = refIndexS_.sequenceMap;

    // Output file stream
    SAVE SAVECLASS(outFileName_);

    stringstream outStream; // Use stringstream instead of string concatenation
    static const int32_t CACHE_SIZE = 1024 * 1024 * 10; // Cache size is 10mb
    outStream.str().reserve(CACHE_SIZE);

    // Record the start and chromosome number of the previous variant
    string preChromosome;
    uint32_t preRefStart = 0;

    // Input file stream
    // Store vcf information
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(vcfFileName_);

    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty line, skip
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }
        
        // comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            if (INFOSTRUCTTMP.line.find("#CHROM") != string::npos) {
                // 1. Head plus chromosome length. (bayestyper)
                outStream << refIndexS_.chrLenTxt;

                // save header
                outStream << INFOSTRUCTTMP.line + "\n";
            }
            // Skip if the vcf contains chromosome length information
            else if (INFOSTRUCTTMP.line.find(",length") == string::npos) {
                outStream << INFOSTRUCTTMP.line + "\n";
            }

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
                map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                    INFOSTRUCTTMP.lineVec
                );

                // If there is only one genotype, skip.
                if (GTVecMapTmp.empty()) {
                    continue;
                }
                
                // The minimum allele frequency and deletion rate were calculated
                tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                    GTVecMapTmp, 
                    INFOSTRUCTTMP.lineVec.size() - 9
                );

                // Did not pass the threshold, straight to the next loop
                if (MAFTMP < MAF_ || MISSRATETMP > MISSRATE_) {
                    continue;
                }
            }
        }


        // Check if there is corresponding chromosome information in the submitted genome
        if (seqMap.find(INFOSTRUCTTMP.CHROM) == seqMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: '"
                << INFOSTRUCTTMP.CHROM 
                << "' is not present in '" << refFileName_ << "'" 
                << endl;
            exit(1);
        }

        // 2. Check if qrySeq contains any '<'. If it does, skip it. (<INS>;<DUP>)
        if (INFOSTRUCTTMP.lineVec[4].find("<") != string::npos) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The query sequence contains the '>' symbol, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << "\t" 
                << INFOSTRUCTTMP.lineVec[4] << endl;
            continue;
        }
        

        // 3. The first variation is greater than the read length. (paragraph)
        if (static_cast<uint32_t>(readLen_) > INFOSTRUCTTMP.POS) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: Start of variation is less than read length, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        

        // 4. The sequence must correspond to reference.
        if (seqMap.at(INFOSTRUCTTMP.CHROM).size() < (INFOSTRUCTTMP.END)) {  // Check that the chromosome length is correct
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: The variant end position is greater than the chromosome length -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            exit(1);
        }

        string trueRefSeq = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
        if (trueRefSeq != INFOSTRUCTTMP.REF) {  // If it is different from the sequence in the reference genome, it is replaced with the sequence in the reference genome
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: sequence difference between refgenome and vcf, replace by refgenome sequence -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;

            INFOSTRUCTTMP.lineVec[3] = trueRefSeq;
            INFOSTRUCTTMP.REF = trueRefSeq;
        }


        // 5. Replace the 'END=' in the eighth column of the VCF. (paragrpah)
        // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end)))
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


        // 6. Graphtyper2 will report an error when it encounters SVLEN=1,2.
        INFOSTRUCTTMP.lineVec[7] = regex_replace(INFOSTRUCTTMP.lineVec[7], svlenReg, "");


        // 7. The first base of refSeq is the same as the first base of qrySeq. (paragraph)
        for (size_t i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++) {
            string qrySeq = INFOSTRUCTTMP.ALTVec[i];

            // When the first base is different, add the preceding base of refSeq to the sequence of qrySeq
            if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1)) {
                // Move the coordinates one step backward
                INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
                INFOSTRUCTTMP.lineVec[1] = to_string(INFOSTRUCTTMP.POS);
                INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
                INFOSTRUCTTMP.REF = seqMap.at(INFOSTRUCTTMP.CHROM).substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // Reextract sequence information
                INFOSTRUCTTMP.lineVec[3] = INFOSTRUCTTMP.REF; // Assigns a value to the Vector

                // Add 'refSeq[0]' to all sequences in 'INFOSTRUCTTMP.ALTVec'
                for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++) {
                    INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                }

                qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq Reassigns a value
            }

            // 8. Check whether ref and qry sequences are the same and the same sites are skipped
            if (qrySeq == INFOSTRUCTTMP.REF) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Sequence same in REF and ALT, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
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
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                INFOSTRUCTTMP.line.clear();
                break;
            }
        }
        // Replace the qry sequence
        INFOSTRUCTTMP.lineVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


        // 10. Check whether the qry after replacement is the same. If yes, skip this site.
        // bayestyper Requirements (A ATG,TGT,ATG)
        // Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed.
        for (const auto& it : INFOSTRUCTTMP.ALTVec) {
            if (count(INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end(), it) > 1) {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: Allelic repeat, skip this site -> "
                    << INFOSTRUCTTMP.CHROM << "\t" 
                    << INFOSTRUCTTMP.POS << endl;
                INFOSTRUCTTMP.line.clear();
                break;
            }
        }


        // 11. Check for mutations that duplicate positions
        // If it is a new chromosome, reset the starting position to zero
        if (INFOSTRUCTTMP.CHROM != preChromosome) {
            preChromosome = INFOSTRUCTTMP.CHROM;
            preRefStart = 0;
        }
        if (INFOSTRUCTTMP.POS == preRefStart) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: multiple variants observed, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // Find the location of gt in the FORMAT field
        vector<string> formatVec = split(INFOSTRUCTTMP.FORMAT, ":");
        vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
        uint32_t maxGT = 0; // Record the largest GT, see if there are more sequences than qry, and skip this site if there are more
        // As long as GT field
        INFOSTRUCTTMP.lineVec[8] = "GT";

        if (gtItera != formatVec.end()) {  // The FORMAT contains GT
            // GT index
            uint32_t gtIndex = distance(formatVec.begin(), gtItera);
            
            // Loop over the GT column
            for (size_t i = 9; i < INFOSTRUCTTMP.lineVec.size(); i++) {
                // Find the genotype
                string gt = split(INFOSTRUCTTMP.lineVec[i], ":")[gtIndex];

                // 12. Combine the genotypes. Convert to.|.. (PanGenie)
                if (gt == ".") {
                    gt = ".|.";
                } else if (gt == "0") {
                    gt = "0|0";
                }
                

                // 13. Change the/in genotype to |. (PanGenie)
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

                // 14. Retain only the diploid variation in the genotype
                if (gtVec.size() == 1) {
                    gt = gtVec[0] + "|0";
                } else if (gtVec.size() > 2) {
                    gt = gtVec[0] + "|" + gtVec[1];
                }
                
                INFOSTRUCTTMP.lineVec[i] = gt;
            }
        }
        else {  // Skip the site
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: GT not in FORMAT column, skip this site -> " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        // 15. Check if GT has more sequences than qry
        // pangenie: VariantReader::VariantReader: invalid genotype in VCF.
        if (INFOSTRUCTTMP.ALTVec.size() < maxGT) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: there are more GTs than qry sequences: " 
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }
        
        // Reassign preRefStart
        preRefStart = INFOSTRUCTTMP.POS;

        // Check whether the information is cleared. If it is cleared, skip it
        if (INFOSTRUCTTMP.line.empty()) {
            continue;
        }

        // Add the replaced string to outStream
        outStream << join(INFOSTRUCTTMP.lineVec, "\t") + "\n";

        if (outStream.tellp() >= CACHE_SIZE) {  // Cache size is 10mb
            string outTxt = outStream.str();
            SAVECLASS.save(outTxt);
            // Clear stringstream
            outStream.str(string());
            outStream.clear();
        }
    }

    if (outStream.tellp() > 0) {  //Write for the last time
        string outTxt = outStream.str();
        SAVECLASS.save(outTxt);
        // Clear stringstream
        outStream.str(string());
        outStream.clear();
    }
}