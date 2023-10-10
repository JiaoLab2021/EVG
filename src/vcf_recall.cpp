// g++ vcf_recall.cpp -o vcf_recall -lz -O2

#include "../include/vcf_recall.hpp"

using namespace std;

int main_recall(int argc, char** argv)
{
    string trueVCFFileName;
    string evaluateFileName;

    string model = "recall";
    string rocKey = "";
    // roc calculation rule
    string rocType = "genotype";

    // Input parameter
    int c;
    while (true)
    {
        static const struct option long_options[] = {
            {"vcf1", required_argument, 0, 't'},
            {"vcf2", required_argument, 0, 'e'},         
            {"model", required_argument, 0, 'm'},

            {"score-field", required_argument, 0, 's'},
            {"roc", required_argument, 0, 'r'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "t:e:m:s:r:h", long_options, &option_index);

        if (c == -1)
            break;

        switch (c)
        {
        case 't':
            trueVCFFileName = optarg;
            break;
        case 'e':
            evaluateFileName = optarg;
            break;
        case 'm':
            model = optarg;
            break;
        case 's':
            rocKey = optarg;
            break;
        case 'r':
            rocType = optarg;
            break;
        case 'h':
        case '?':
            help_recall(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // Check whether the parameters are correct
    if (trueVCFFileName.empty() || evaluateFileName.empty() || model.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_recall(argv);
        return 1;
    }

    // Check that rocType is correct
    if (rocType != "recall" && rocType != "genotype") {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error: -r\n";
        help_recall(argv);
        return 1;
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Running ...\n";

    // init
    VCFRecall VCFRecallClass(
        trueVCFFileName, 
        evaluateFileName, 
        rocKey, 
        rocType
    );
    // index
    VCFRecallClass.build_true_index();


    if (model == "recall") {
        VCFRecallClass.evulate_recall();
    } else if (model == "genotype") {
        VCFRecallClass.evulate_gt();
    } else if (model == "gramtools") {
        VCFRecallClass.gramtools_convert();
        VCFRecallClass.evulate_gt();
    } else {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "The model is incorrect. [recall/genotype/gramtools]" << endl;
        exit(1);
    }
    
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Done ...\n";

    return 0;
}

// Help document
void help_recall(char** argv)
{
  cerr << "usage: " << argv[0] << " " << argv[1] << " -t -e [options]" << endl
       << "evaluate the results of genome graph software." << endl
       << endl
       << "input files (Note: vcf files must be sorted):" << endl
       << "    -t, --vcf1          FILE       vcf file containing truth values" << endl
       << "    -e, --vcf2          FILE       output vcf file by genomes graph software" << endl
       << endl
	   << "algorithm:" << endl
       << "    -m, --model         STRING     the mode in which the software runs (recall/genotype/gramtools) [recall]" << endl
       << endl 
       << "ROC arguments:" << endl
       << "    -s, --score-field   STRING     the name of the VCF FORMAT/INFO field to use as the ROC score (FORMAT.<name>/INFO.<name>) [FORMAT.DP]" << endl
       << "    -r, --roc           STRING     roc calculation rules (recall/genotype) [genotype]" << endl
       << endl
       << "    -h, --help                     print this help document" << endl;
}


/**
 * init
 *
 * @param trueVCFFileName    true set
 * @param evaluateFileName   evaluation set
 * @param rocKey             Keywords used to extract roc score
 * 
**/
VCFRecall::VCFRecall(
    const string& trueVCFFileName, 
    const string& evaluateFileName, 
    const string& rocKey, 
    const string& rocType
) : trueVCFFileName_(trueVCFFileName), evaluateFileName_(evaluateFileName), rocKey_(rocKey), rocType_(rocType)
{
    // Variable length range
    lengthVec_.push_back("-999999/-10000");
    lengthVec_.push_back("-9999/-5000");
    lengthVec_.push_back("-4999/-2500");
    lengthVec_.push_back("-2499/-1000");
    lengthVec_.push_back("-999/-500");
    lengthVec_.push_back("-499/-400");
    lengthVec_.push_back("-399/-300");
    lengthVec_.push_back("-299/-200");
    lengthVec_.push_back("-199/-100");
    lengthVec_.push_back("-99/-50");
    lengthVec_.push_back("-49/-1");
    lengthVec_.push_back("0/0");
    lengthVec_.push_back("1/49");
    lengthVec_.push_back("50/99");
    lengthVec_.push_back("100/199");
    lengthVec_.push_back("200/299");
    lengthVec_.push_back("300/399");
    lengthVec_.push_back("400/499");
    lengthVec_.push_back("500/999");
    lengthVec_.push_back("1000/2499");
    lengthVec_.push_back("2500/4999");
    lengthVec_.push_back("5000/9999");
    lengthVec_.push_back("10000/999999");

    // ROC -> save the number of DP or GQ
    rocColNum_ = 0;
    if (!rocKey_.empty())
    {
        rocKeyVec_ = split(rocKey_, "."); // vector<colname, key>

        if (rocKeyVec_[0] == "INFO")
        {
            rocColNum_ = 7;
        }
        else if (rocKeyVec_[0] == "FORMAT")
        {
            rocColNum_ = 8;
        }
    }
}


/**
 * build the index of true VCF file
 * 
 * @return void
**/
void VCFRecall::build_true_index()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index ...\n";

    // Check if the vcf file is sorted
    check_vcf_sort(trueVCFFileName_);

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(trueVCFFileName_);
    
    // If not traversed, continue
    while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // empty or commant line, skip
        if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );

        string gt = join(gtVec, "/");

        // Help document
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  // If it is (0/0,.) Format is skipped, or returns an empty list without building an index
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The genotype of the variant is empty, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // Gets the length of the ref and haplotype
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // Gets the length of the variation
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // Save the length of all variations
        trueVCFStrcuture_.allLengthList.push_back(svLength);

        // Variation information
        trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);

        // Initialize the hash table first
        if (trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM) == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM];
        }
        trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS] = make_tuple(refLen, qryLenVec, INFOSTRUCTTMP.line, svLength, gtVec);
    }
}


/**
 * Convert the results of gramtools
 * 
 * @return void
**/
void VCFRecall::gramtools_convert()
{
    string outInformations;

    // check whether the file is sorted
    check_vcf_sort(evaluateFileName_);

    // output file stream
    vector<string> prefixVec = split(evaluateFileName_, "/");
    string outFileName = "convert." + prefixVec.back();
    SAVE SaveClass(outFileName);

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);
    // If not traversed, continue
    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            SaveClass.save(INFOSTRUCTTMP.line);
        }

        // Extract genotype and locus coverage information from informationsVector
        vector<string> gtVector = split(INFOSTRUCTTMP.lineVec.back(), ":");
        vector<string> formatVector = split(INFOSTRUCTTMP.FORMAT, ":");
        string gt;
        if (gtVector.size() > 0) {
            gt = gtVector[0];
        } else {
            gt = "0";
        }

        // Find the location of COV in the FORMAT field
        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
        int covIndex = 0;
        if (covItera != formatVector.end()) { // The FORMAT contains COV
            covIndex = distance(formatVector.begin(), covItera);
        } else {  // If not, just switch and continue the next cycle
            outInformations += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // The contents of the FILTER field are added directly to outInformations
        for (size_t i = 0; i < INFOSTRUCTTMP.lineVec.size()-1; i++) {
            if (i == 6) {  // The FILTER column is changed to PASS, and the source file is.
                outInformations += "PASS\t";
            } else {
                outInformations += INFOSTRUCTTMP.lineVec[i] + "\t";
            }
        }

        string siteCoverage = gtVector[covIndex];
        
        // If gt=0 means that there is no variation at this site, then gt is changed to 0/0
        if (gt == "0") {
            outInformations += "\t0/0";
        } else {  // When gt is not 0
            // If there is a 0 field in the site coverage, it represents a mutation site that is a pure sum
            if (siteCoverage.find("0,") < siteCoverage.length()) {
                outInformations += "\t1/1";
            } else {  // If the site coverage does not have a 0, field, it indicates that the site is a heterozygous mutation site
                outInformations += "\t0/1";
            }
        }

        // Add the last field of informations
        for (size_t i = 1; i < gtVector.size(); i++) {
            outInformations += ":" + gtVector[i];
        }

        outInformations += "\n";

        // The file is written every 10Mb
        if (outInformations.size() > bufferSize_) {
            SaveClass.save(outInformations);

            // The file is written every 10Mb
            outInformations.clear();
        }
    }

    // Last write to the file
    SaveClass.save(outInformations);
    
    evaluateFileName_ = outFileName;
}


/**
 * genotype evaluation function
 * 
 * @return void
**/
void VCFRecall::evulate_gt()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Evaluating the genotype results ..." << endl;

    // check whether the file is sorted
    check_vcf_sort(evaluateFileName_);
    
    // Save the correct results for typing
    stringstream truetxtStream; // Use stringstream instead of string concatenation
    truetxtStream.str().reserve(bufferSize_);
    SAVE trueFile("genotype.true.vcf.gz");

    // Saves the result of a typing error
    stringstream genotypeMisTxtStream; // Use stringstream instead of string concatenation
    genotypeMisTxtStream.str().reserve(bufferSize_);
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // Save the result of miscall
    stringstream misCallTxtStream; // Use stringstream instead of string concatenation
    misCallTxtStream.str().reserve(bufferSize_);
    SAVE misCallFile("miscall.vcf.gz");

    // No variation found
    stringstream failCallTxtStream; // Use stringstream instead of string concatenation
    failCallTxtStream.str().reserve(bufferSize_);
    SAVE failCallFile("failcall.vcf.gz");

    // No variation found
    vector<int64_t> genotypeLenVec;
    vector<int64_t> misgenotypeLenVec;
    vector<int64_t> miscallLenVec;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // Skip comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }
        
        // According to the FILTER field
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // If not (PASS), skip directly
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // Get genotype information
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // Filter by genotype
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  // If it is (0/0,.) Format is skipped, or returns an empty list without building an index
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The genotype of the variant is empty, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        /* --------------------------------- roc ------------------------------------ */
        float rocNum = 0.0;
        rocNum = get_roc(INFOSTRUCTTMP);
        /* --------------------------------- roc ------------------------------------ */

        // Gets the length of the ref and haplotype
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // Gets the length of the variation
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // Add length to rocAllMap
        rocCallMap_[rocNum].push_back(svLength);

        // Genotyping was performed
        // Check for the presence of chromosomes. No is a miscall
        auto iter1 = trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM);
        if (iter1 == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: " << INFOSTRUCTTMP.CHROM << " not in true set.\n";
            miscallLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // Record the result of the genotype
        uint32_t genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype
        uint32_t trueRefLen;
        string trueVcfInfo;
        int32_t trueSvLen;

        // Find the true refLen in the map of the corresponding chromosome qryLenVec
        auto iter2 = iter1->second.find(INFOSTRUCTTMP.POS);
        if (iter2 != iter1->second.end()) {  // The software finds it and then determines whether the classification is correct
            // true mutation information
            trueRefLen = get<0>(iter2->second);
            const vector<uint32_t>& trueQryLenVec = get<1>(iter2->second);
            trueVcfInfo = get<2>(iter2->second);
            trueSvLen = get<3>(iter2->second);
            vector<int> trueGtVec = get<4>(iter2->second);

            // If the length is 0, an error is reported and the code exits.
            if (trueQryLenVec.size() == 0) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: trueQryLenVec.size() == 0 -> " 
                    << INFOSTRUCTTMP.CHROM << " " 
                    << INFOSTRUCTTMP.POS << endl;
                exit(1);
            }

            // Traverse alleles
            int64_t trueQryLenVecIdx = -1;  // Record used haplotypes for truth values to prevent duplicate evaluations
            for (size_t i = 0; i < qryLenVec.size(); i++) {
                uint32_t qryLen = qryLenVec[i];
                int callGt = gtVec[i];

                for (size_t j = 0; j < trueQryLenVec.size(); j++) {
                    // Used haplotypes are skipped
                    if (static_cast<int64_t>(j) == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[j];
                    int trueGt = trueGtVec[j];

                    // If one of callGt and trueGt is 0, but the other is not, the next loop will be performed to prevent SNP from making errors in judgment
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // deficiency
                    if (refLen >= 50 && qryLen < 50) {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // ins
                        if ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // Replace
                        if (((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {  // snp
                        if (((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {  // del
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen <= 2 && qryLen >= 3) {  // ins
                        if ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25 && 
                            (std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    }
                }
                stop:; // If found, exit the nested loop. End the cycle of this haplotype and proceed to the next haplotype
            }
        }

        // Determine the result of typing
        if (genotypeTrueNum == 0) {  // The mutations we found are not in the true concentration 
            miscallLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        } else {
            // Correct typing
            if (genotypeTrueNum >= gtVec.size()) {
                genotypeLenVec.push_back(trueSvLen);
                truetxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // Add length to rocTrueMap
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // Typing error
                misgenotypeLenVec.push_back(trueSvLen);
                genotypeMisTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // The roc is calculated using the recall, or only the genotype if the roc is not enabled
                if (rocType_ == "recall") {
                    // Add length to rocTrueMap
                    rocRecallMap_[rocNum].push_back(svLength);
                }
                
            }
            // Delete variations that have already been used
            trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
        }

        // save result
        if (truetxtStream.tellp() >= bufferSize_ || 
            genotypeMisTxtStream.tellp() >= bufferSize_ || 
            misCallTxtStream.tellp() >= bufferSize_) // It is written every 10Mb
        {
            // Correct genotyping
            string truetxt = truetxtStream.str();
            trueFile.save(truetxt);
            // Clear stringstream
            truetxtStream.str(string());
            truetxtStream.clear();

            // genotyping error
            string genotypeMisTxt = genotypeMisTxtStream.str();
            genotypeMisFile.save(genotypeMisTxt);
            // Clear stringstream
            genotypeMisTxtStream.str(string());
            genotypeMisTxtStream.clear();

            // Variation that is not in true concentration
            string misCallTxt = misCallTxtStream.str();
            misCallFile.save(misCallTxt);
            // Clear stringstream
            misCallTxtStream.str(string());
            misCallTxtStream.clear();
        }
    }

    // last save
    if (truetxtStream.tellp() > 0 || 
        genotypeMisTxtStream.tellp() > 0 || 
        misCallTxtStream.tellp() > 0) // It is written every 10Mb
    {
        // Correct typing
        string truetxt = truetxtStream.str();
        trueFile.save(truetxt);
        // Clear stringstream
        truetxtStream.str(string());
        truetxtStream.clear();

        // genotyping error
        string genotypeMisTxt = genotypeMisTxtStream.str();
        genotypeMisFile.save(genotypeMisTxt);
        // Clear stringstream
        genotypeMisTxtStream.str(string());
        genotypeMisTxtStream.clear();

        // Variation that is not in true concentration
        string misCallTxt = misCallTxtStream.str();
        misCallFile.save(misCallTxt);
        // Clear stringstream
        misCallTxtStream.str(string());
        misCallTxtStream.clear();
    }

    // True variation not found
    for (auto it1 : trueVCFStrcuture_.chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (auto it2 : it1.second) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            failCallTxtStream << get<2>(it2.second) + "\n";

            // save the result
            if (failCallTxtStream.tellp() >= bufferSize_) {  // It is written every 10Mb
                string failCallTxt = failCallTxtStream.str();
                failCallFile.save(failCallTxt);
                // Clear stringstream
                failCallTxtStream.str(string());
                failCallTxtStream.clear();
            }
        }
    }
    // last save
    if (failCallTxtStream.tellp() > 0) {  // It is written every 10Mb
        string failCallTxt = failCallTxtStream.str();
        failCallFile.save(failCallTxt);
        // Clear stringstream
        failCallTxtStream.str(string());
        failCallTxtStream.clear();
    }

    int64_t sv_genotype_recall = genotypeLenVec.size();
    int64_t sv_misgenotype_recall = misgenotypeLenVec.size();
    int64_t sv_mis_call = miscallLenVec.size();
    int64_t sv_all = trueVCFStrcuture_.allLengthList.size();

    // save result
    ofstream outFile;
    outFile.open("vcf_evulate.out", ios::app);

    outFile << "snp+indel+sv:\n"
            << "genotype_recall:" << sv_genotype_recall 
            << "\nmisgenotype_call:" << sv_misgenotype_recall 
            << "\nrecall:" << sv_genotype_recall + sv_misgenotype_recall 
            << "\nmis_call:" << sv_mis_call 
            << "\ncall:" << sv_genotype_recall + sv_misgenotype_recall + sv_mis_call 
            << "\nfail_call:" << sv_all - sv_genotype_recall - sv_misgenotype_recall 
            << "\nall:" << sv_all 
            << endl
            << endl;

    vector<int64_t> all_length_count = count_num(trueVCFStrcuture_.allLengthList);
    vector<int64_t> mis_call_length_count = count_num(miscallLenVec);
    vector<int64_t> genotype_length_count = count_num(genotypeLenVec);
    vector<int64_t> misgenotype_length_count = count_num(misgenotypeLenVec);
    
    outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

    // True variation not found
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        outFile << lengthVec_[i] << ": " 
                << genotype_length_count[i] << "/" 
                << misgenotype_length_count[i] << "/" 
                << genotype_length_count[i] + misgenotype_length_count[i] << "/" 
                << mis_call_length_count[i] << "/" 
                << genotype_length_count[i] + misgenotype_length_count[i] + mis_call_length_count[i] << "/" 
                << all_length_count[i] - genotype_length_count[i] - misgenotype_length_count[i] << "/" 
                << all_length_count[i] << endl;
        
        if ( 10 <= i && i <= 12)
        {
            sv_genotype_recall -= genotype_length_count[i];
            sv_misgenotype_recall -= misgenotype_length_count[i];
            sv_mis_call -= mis_call_length_count[i];
            sv_all -= all_length_count[i];
        }
    }

    outFile << "\nsv:\n"
            << "genotype_recall:" << sv_genotype_recall 
            << "\nmisgenotype_call:" << sv_misgenotype_recall 
            << "\nrecall:" << sv_genotype_recall + sv_misgenotype_recall 
            << "\nmis_call:" << sv_mis_call 
            << "\ncall:" << sv_genotype_recall + sv_misgenotype_recall + sv_mis_call 
            << "\nfail_call:" << sv_all - sv_genotype_recall - sv_misgenotype_recall 
            << "\nall:" << sv_all 
            << endl
            << endl;

    // Free memory
    outFile.close();

    roc_calculate(trueVCFStrcuture_.allLengthList);
}


/**
 * genotype evaluation function
 * 
 * @return void
**/
void VCFRecall::evulate_recall()
{
    cerr << "[" << __func__ << "::" << getTime() << "] " << "Evaluating the genotype results ..." << endl;

    // Check whether the vcf file is sorted
    check_vcf_sort(evaluateFileName_);

    // Save the correct results
    stringstream truetxtStream; // Use stringstream instead of string concatenation
    truetxtStream.str().reserve(bufferSize_);
    SAVE trueFile("recall.true.vcf.gz");

    // Save the correct results for typing
    stringstream trueGtTxtStream; // Use stringstream instead of string concatenation
    trueGtTxtStream.str().reserve(bufferSize_);
    SAVE trueGtFile("genotype.true.vcf.gz");

    // Saves the result of a typing error
    stringstream genotypeMisTxtStream; // Use stringstream instead of string concatenation
    genotypeMisTxtStream.str().reserve(bufferSize_);
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // Save the result of miscall
    stringstream misCallTxtStream; // Use stringstream instead of string concatenation
    misCallTxtStream.str().reserve(bufferSize_);
    SAVE misCallFile("miscall.vcf.gz");

    // Define vector
    vector<int64_t> genotypeLenVec;
    vector<int64_t> callLenVec;
    vector<int64_t> recallLenVec;

    string preChr; // It's used to determine if it's a new chromosome, and if it is, it reconstructs the vector.

    // Binary search method left and right index
    int leftIdxTmp = 0;
    int rightIdxTmp = 0;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // Skip comment line
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        // According to the FILTER field
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // Skip comment line
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // According to the FILTER field
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // Filter by genotype
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  // If it is (0/0,.) Format is skipped, or returns an empty list without building an index
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The genotype of the variant is empty, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }


        /* --------------------------------- roc ------------------------------------ */
        float rocNum = 0.0;
        rocNum = get_roc(INFOSTRUCTTMP);
        /* --------------------------------- roc ------------------------------------ */


        // Gets the length of the ref and haplotype
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // Gets the length of the variation
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // Add length to rocAllMap
        rocCallMap_[rocNum].push_back(svLength);


        // recall
        // Check for the presence of chromosomes. No is a miscall
        if (trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM) == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.CHROM << "' is not in the true set.\n";
            callLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }
    

        if (preChr != INFOSTRUCTTMP.CHROM) {  // Check for the presence of chromosomes. No is a miscall
            // The left and right indexes of binary search are zeroed out to speed up the query
            leftIdxTmp = 0;
            rightIdxTmp = 0;
        }

        // Record binary search index, multiple alleles with the same refStart
        int64_t indexLeft = -1;
        int64_t indexRight = -1;

        uint32_t genotypeTrueNum = 0;  // 0->misCall >0 && <gtVec.size()->recall  >=gtVec.size()->trueGenotype

        // recalled  trueRefStart, trueRefLen, trueSvLen°¢truevcfInfo
        uint32_t trueRefStart;
        uint32_t trueRefLen;
        int32_t trueSvLen;
        string truevcfInfo;

        // Record the used haplotypes in the true set
        int64_t trueQryLenVecIdx = -1;

        // Traverse multiple alleles
        for (size_t i = 0; i < qryLenVec.size(); i++) {
            uint32_t qryLen = qryLenVec[i];
            int callGt = gtVec[i];

            // Binary search method to find the variation index within 200/10bp.
            if (indexLeft == -1 && indexRight == -1) {  // New alleles were searched again
                if (refLen <= 49 && qryLen <= 49) {
                    indexLeft = search_Binary_left(trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS-10, leftIdxTmp);
                    leftIdxTmp = indexLeft;
                    indexRight = search_Binary_right(trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS+10, rightIdxTmp);
                    rightIdxTmp = indexRight;
                } else {
                    indexLeft = search_Binary_left(trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS-200, leftIdxTmp);
                    leftIdxTmp = indexLeft;
                    indexRight = search_Binary_right(trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM], INFOSTRUCTTMP.POS+200, rightIdxTmp);
                    rightIdxTmp = indexRight;
                } if (indexLeft < 0 || indexRight >= static_cast<int64_t>(trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].size())) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: out of index, please check the data or code.\n";
                    exit(1);
                }
            }

            // Traverse the index of the binary search method
            for (int64_t j = indexLeft; j <= indexRight; j++) {
                // Traverse the index of the binary search method
                trueRefStart = trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM][j];

                // First check to see if the mutation has been used, then move on to the next coordinates if it's been deleted
                if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(trueRefStart) == trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                    continue;
                }
                
                trueRefLen = get<0>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                const vector<uint32_t>& trueQryLenVec = get<1>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                truevcfInfo = get<2>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                trueSvLen = get<3>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                vector<int> trueGtVec = get<4>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);

                for (size_t k = 0; k < trueQryLenVec.size(); k++) {  // First check to see if the mutation has been used, then move on to the next coordinates if it's been deleted
                    // Used haplotypes are skipped
                    if (static_cast<int64_t>(k) == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[k];
                    int trueGt = trueGtVec[k];

                    // If one of callGt and trueGt is 0, but the other is not, the next loop will be performed to prevent SNP from making errors in judgment
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // deletion
                    if (refLen >= 50 && qryLen < 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // The haplotype recall and genotype were recorded correctly
                            genotypeTrueNum++;

                            // Multiple alleles use the same starting position
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // ≤Â»Î
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // The haplotype recall and genotype were recorded correctly
                            genotypeTrueNum++;

                            // Multiple alleles use the same starting position
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // Replace
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // The haplotype recall and genotype were recorded correctly
                            genotypeTrueNum++;

                            // Multiple alleles use the same starting position
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {  // snp
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // Records indexes that have used haplotypes
                            genotypeTrueNum++;

                            // Multiple alleles use the same starting position
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {  // indel-del
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // The haplotype recall and genotype were recorded correctly
                            genotypeTrueNum++;

                            // Multiple alleles use the same starting position
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }  else if (refLen <= 2 && qryLen >= 3) {  // indel-ins
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // The haplotype recall and genotype were recorded correctly
                            genotypeTrueNum++;

                            // The haplotype recall and genotype were recorded correctly
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // Records indexes that have used haplotypes
                            genotypeTrueNum++;

                            // The haplotype recall and genotype were recorded correctly
                            indexLeft = j;
                            indexRight = j;

                            // Records indexes that have used haplotypes
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }
                }
            }
            stop:; // If found, exit the nested loop. End the cycle of this haplotype and proceed to the next haplotype
        }

        // Determine the results of the search and add 
        if (genotypeTrueNum > 0) {  // Determine the results of the search and add
            // Software call vcf-vector
            callLenVec.push_back(trueSvLen);

            truetxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
            recallLenVec.push_back(trueSvLen);

            // Correct typing
            if (genotypeTrueNum >= gtVec.size()) {  // If the value is greater than or equal to gtVec.size(), the classification is correct
                genotypeLenVec.push_back(trueSvLen);
                trueGtTxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                // Add length to rocTrueMap
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // Typing error
                genotypeMisTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                // The roc is calculated using the recall, or only the genotype if the roc is not enabled
                if (rocType_ == "recall")
                {
                    // Add length to rocTrueMap
                    rocRecallMap_[rocNum].push_back(svLength);
                }
            }

            // Add length to rocTrueMap
            if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(INFOSTRUCTTMP.POS) != trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
            }
        } else {  // If the above loop does not find the variation in the true set, it is added with the length found by the software itself
            callLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
        }

        // save result
        if (truetxtStream.tellp() >= bufferSize_ || 
            trueGtTxtStream.tellp() >= bufferSize_ || 
            genotypeMisTxtStream.tellp() >= bufferSize_ || 
            misCallTxtStream.tellp() >= bufferSize_) // It is written every 10Mb
        {
            // Correct call
            string truetxt = truetxtStream.str();
            trueFile.save(truetxt);
            // Clear stringstream
            truetxtStream.str(string());
            truetxtStream.clear();

            // Correct genotyping
            string trueGtTxt = trueGtTxtStream.str();
            trueGtFile.save(trueGtTxt);
            // Clear stringstream
            trueGtTxtStream.str(string());
            trueGtTxtStream.clear();

            // Genotyping error
            string genotypeMisTxt = genotypeMisTxtStream.str();
            genotypeMisFile.save(genotypeMisTxt);
            // Clear stringstream
            genotypeMisTxtStream.str(string());
            genotypeMisTxtStream.clear();

            // Variation that is not in true concentration
            string misCallTxt = misCallTxtStream.str();
            misCallFile.save(misCallTxt);
            // Clear stringstream
            misCallTxtStream.str(string());
            misCallTxtStream.clear();
        }
    }

    // last save
    if (truetxtStream.tellp() > 0 || 
        trueGtTxtStream.tellp() > 0 || 
        genotypeMisTxtStream.tellp() > 0 || 
        misCallTxtStream.tellp() > 0) // It is written every 10Mb
    {
        // Correct call
        string truetxt = truetxtStream.str();
        trueFile.save(truetxt);
        // Clear stringstream
        truetxtStream.str(string());
        truetxtStream.clear();

        // Correct genotyping
        string trueGtTxt = trueGtTxtStream.str();
        trueGtFile.save(trueGtTxt);
        // Clear stringstream
        trueGtTxtStream.str(string());
        trueGtTxtStream.clear();

        // Genotyping error
        string genotypeMisTxt = genotypeMisTxtStream.str();
        genotypeMisFile.save(genotypeMisTxt);
        // Clear stringstream
        genotypeMisTxtStream.str(string());
        genotypeMisTxtStream.clear();

        // Variation that is not in true concentration
        string misCallTxt = misCallTxtStream.str();
        misCallFile.save(misCallTxt);
        // Clear stringstream
        misCallTxtStream.str(string());
        misCallTxtStream.clear();
    }

    // Save no variation found, false negative
    saveFailCall(
        trueVCFStrcuture_.chrStartLenInfoGtTupMap, 
        "failcall.vcf.gz"
    );

    int64_t all_num = trueVCFStrcuture_.allLengthList.size();
    int64_t call_num = callLenVec.size();
    int64_t recall_num = recallLenVec.size();

    int64_t genotype_recall_number = genotypeLenVec.size();

    // save result
    ofstream outFile;
    outFile.open("vcf_evulate.out", ios::app);

    outFile << "snp+indel+sv:\n"
            << "genotype_recall:" << genotype_recall_number 
            << "\nmisgenotype_call:" << recall_num - genotype_recall_number 
            << "\nrecall:" << recall_num 
            << "\nmis_call:" << call_num - recall_num 
            << "\ncall:" << call_num 
            << "\nfail_call:" << all_num - recall_num 
            << "\nall:" << all_num 
            << endl 
            << endl;

    vector<int64_t> allNumVec = count_num(trueVCFStrcuture_.allLengthList);
    vector<int64_t> callNumVec = count_num(callLenVec);
    vector<int64_t> recallNumVec = count_num(recallLenVec);
    vector<int64_t> genotypeRecallNumVec = count_num(genotypeLenVec);
    
    outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

    // Save no variation found, false negative
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        outFile << lengthVec_[i] << ": " 
                << genotypeRecallNumVec[i] << "/" 
                << recallNumVec[i] - genotypeRecallNumVec[i] << "/" 
                << recallNumVec[i] << "/" 
                << callNumVec[i] - recallNumVec[i] << "/"
                << callNumVec[i] << "/"
                << allNumVec[i] - recallNumVec[i] << "/"
                << allNumVec[i] << endl;
        
        if (10 <= i && i <= 12) {
            all_num -= allNumVec[i];
            call_num -= callNumVec[i];
            recall_num -= recallNumVec[i];
            genotype_recall_number -= genotypeRecallNumVec[i];
        }    
    }

    outFile << "\nsv:\n"
            << "genotype_recall:" << genotype_recall_number 
            << "\nmisgenotype_call:" << recall_num - genotype_recall_number 
            << "\nrecall:" << recall_num 
            << "\nmis_call:" << call_num - recall_num 
            << "\ncall:" << call_num 
            << "\nfail_call:" << all_num - recall_num 
            << "\nall:" << all_num 
            << endl 
            << endl;

    // Length result output
    outFile.close();

    // Free memory
    roc_calculate(trueVCFStrcuture_.allLengthList);
}


/**
 * Get a list of loci genotypes.
 *
 * @param lineVec          lineVec
 * @param sampleIdx        Index of the sample genotype, with the default value 0 representing the last column
 * 
 * 
 * @return gtVec           vector <int>
**/
vector<int> VCFRecall::get_gt(
    const vector<string> & lineVec, 
    int sampleIdx
) {
    vector <int> gtVec;  // Locus typing vector

    // FORMAT character split
    vector<string> formatVec = split(lineVec[8], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // Gets the index location of GT

    if (gtIndex == formatVec.size()) {  // Gets the index location of GT
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No [GT] information in FORMAT, replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {  // If it exists, save it
        string gt;  // Store the genotype field

        if (sampleIdx == 0) {  // The number of columns not specified is the last column
            gt = split(lineVec.back(), ":")[gtIndex];  // gt field
        } else {
            gt = split(lineVec[sampleIdx], ":")[gtIndex];  // gt field
        }
        
        string splitStr;  // gt field
        if (gt.find("/") != string::npos) {  // Delimiter in gt
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {  // Determine that '|' is the separator
            splitStr = "|";
        } else {  // If the delimiter is not known, use 0|0 instead
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: Genotype does not contain '\\' or '|', replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
            gtVec = {0, 0};
            return gtVec;
        }

        // Once you find gt, split it by splitStr and loop
        auto splitResult = split(gt, splitStr);
        for (auto it : splitResult) {
            // If it is '.', skip this site
            if (it == ".") {
                gtVec = {0, 0};
                return gtVec;
            }
            // Prevent error in stof(it) conversion
            try {
                gtVec.push_back(stoi(it));  // Add to vector
            } catch (const std::exception& e) {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: stof(it) conversion failed. Replaced with {0,0} -> " << it << endl;
                gtVec = {0, 0};
                return gtVec;
            }
        }
    }

    return gtVec;
}


/**
 * Get a list of loci genotypes.
 *
 * @param INFOSTRUCTTMP    line information
 * 
 * 
 * @return rocNum
**/
float VCFRecall::get_roc(
    const VCFINFOSTRUCT& INFOSTRUCTTMP
)
{
    if (rocKeyVec_.size() < 2) {
        return 0.0;
    }

    if (rocColNum_ == 8) {  // FORMAT field
        vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec.back(), ":");

        vector<string> rocInfVec = split(INFOSTRUCTTMP.FORMAT, ":");
        // rocNum subscript
        vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec_[1]);

        if (rocItera == rocInfVec.end()) {  // rocNum subscript
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: " << rocKeyVec_[1] << " not in " << rocKeyVec_[0] << " column -> " 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        return stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
    } else {  // INFO field
        smatch patternResult; // INFO field
        regex pattern(rocKeyVec_[1] + "=(\\d+)");

        if (!regex_search(INFOSTRUCTTMP.INFO, patternResult, pattern)) {  // Check whether the value of rocKeyVec_[1] is matched
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: " << rocKeyVec_[1] << " not in " << rocKeyVec_[0] << " column -> " 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        string rocNumString = patternResult[1];

        return stof(rocNumString);
    }
}


/**
 * Get the length information of the variation
 *
 * @param refLen            ref length
 * @param qryLenVec         qry Length list
 * @param gtVec             qry Length list
 * 
 * 
 * @return int32_t          svLength
**/
int32_t VCFRecall::sv_length_select(
    const uint32_t & refLen, 
    const vector<uint32_t> & qryLenVec, 
    const vector<int> & gtVec
) {
    int32_t svLength = 0;

    for (size_t i = 0; i < gtVec.size(); i++) {
        if (gtVec[i] == 0) {  // Genotype list
            continue;
        } else {
            svLength = qryLenVec[i] - refLen;
        }
    }
    
    return svLength;
}


/**
 * If genotype is 0, skip
 *
 * @param lengthVec_         Divided interval
 * @param length_list        Length list
 * 
 * 
 * @return vector<int64_t>  The length of each interval
**/
vector<int64_t> VCFRecall::count_num(
    vector<int64_t> length_list
)
{
    vector<int64_t> out_length;
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        int64_t x1 = std::stoll(split(lengthVec_[i], "/")[0]);
        int64_t x2 = std::stoll(split(lengthVec_[i], "/")[1]);

        vector<int64_t>::size_type result = count_if(length_list.begin(), length_list.end(), f_mod(x1, x2));

        out_length.push_back(result);
    }
    return out_length;
}


/**
 * The length of each interval.
 *
 * @param svType                         Variation type
 * @param refSeq                         ref column information
 * @param qrySeqs                        qry Column information
 * @param gtVec                          qry Column information
 * @param lenType                        Take only the length corresponding to the haplotype or all the lengths (hap/all)
 * 
 * 
 * @return tuple<uint32_t, vector<uint32_t> >      tuple<refLen, vector<qryLen> >
**/
tuple<uint32_t, vector<uint32_t> > VCFRecall::get_hap_len(
    const string & svType, 
    const string & refSeq, 
    const string & qrySeqs, 
    const vector<int> & gtVec, 
    const string & lenType
) {
    // Take only the length corresponding to the haplotype or all the lengths (hap/all)
    if (lenType != "hap" && lenType != "all") {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "Error: lenType -> " 
            << lenType << endl;
        exit(1);
    }
    
    uint32_t refLen = refSeq.size();  // Check that the mode is correct and exit the code incorrectly
    vector<string> qrySeqVec = split(qrySeqs, ",");  // ref sequence length
    
    // Check whether the index is out of bounds. If only hap length is needed, check again
    if (lenType == "hap") {
        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
        if (static_cast<uint32_t>(maxGtNum) > qrySeqVec.size()) {  // First check whether the array is out of bounds
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: number of genotyping and ALT sequences do not match -> " 
                << qrySeqs << endl;
            exit(1);
        }
    }

    // First check whether the array is out of bounds
    vector<uint32_t> seqLenVec;
    seqLenVec.push_back(refLen);

    // Construct the length index common to ref and qry
    for (size_t i = 0; i < qrySeqVec.size(); i++) {
        string qrySeq = qrySeqVec[i];

        // Temporary storage length
        uint32_t refLenTmp = refLen;
        uint32_t qryLenTmp = qrySeq.size(); 

        if (qrySeq.find(">") != string::npos) {  // GraphTyper2 results
            // Use regular expressions to extract length information
            std::regex reg("SVSIZE=(\\d+)");
            std::smatch match;
            // Check whether there is any length information in the qry sequence. If there is no length information, skip this site
            if (std::regex_search(qrySeq, match, reg)) {  // <INS:SVSIZE=97:BREAKPOINT1>
                // <INS:SVSIZE=90:BREAKPOINT1>
                if (qrySeq.find("<INS") != string::npos && qrySeq.find("<DEL") == string::npos) {  // The character segment contains inserts
                    refLenTmp = refLen;
                    qryLenTmp = std::stoul(match[1].str());;
                } else if (qrySeq.find("<INS") == string::npos && qrySeq.find("<DEL") != string::npos) {  // The character segment contains missing characters: <DEL:SVSIZE=1720:BREAKPOINT>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLen;
                } else if (qrySeq.find("<DUP") != string::npos) {  // Character segment contains duplicate fields (GraphTyper2) <DUP:SVSIZE=2806:BREAKPOINT1>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLenTmp * 2;
                }
            } else {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No length information in column 4 -> " << qrySeq << endl;
                refLenTmp = refLen;
                qryLenTmp = qrySeq.size();
            }

            seqLenVec[0] = refLenTmp; // Reset the length of the ref
        }

        // The result of BayesTyper was judged to be duplication, ref_len was 1-2, and qry_seq was very long
        if (svType.find("Duplication") != string::npos && refLenTmp <= 2) {  // BayesTyper will change duplication into insertion, so recalculate the length
            refLenTmp = qryLenTmp;
            qryLenTmp *= 2;
        }

        // Add the length of qry
        if (seqLenVec.size() == (i + 1)) {  // Judging that the haplotype is indeed in its place,
            seqLenVec.push_back(qryLenTmp);  // Add the length of qry
        } else {  // Index and list length do not match The Times error
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: Incorrect index for haplotype length -> " << qrySeqs << endl;
            exit(1);
        }
    }

    // Index and list length do not match The Times error
    vector<uint32_t> qryLenVec;
    if (lenType == "hap") {  // The index and list length do not match the qry length corresponding to the genotype found by Times Error a
        if (gtVec.size() == 1) {  // Take only the sequence length of the haplotype
            if (gtVec[0] == 0) {  // If genotype 0, skip this locus
                // Return '0/0'
                vector<uint32_t> qryLenVecTmp(seqLenVec[0], qrySeqVec.size());
                return make_tuple(seqLenVec[0], qryLenVecTmp);
            } else {
                qryLenVec.push_back(seqLenVec[gtVec[0]]);
                qryLenVec.push_back(seqLenVec[gtVec[0]]);
            }
        } else {
            for (auto gtTmp : gtVec) {
                qryLenVec.push_back(seqLenVec[gtTmp]);
            }
        }
    } else {  // Take all the lengths
        for (size_t i = 1; i < seqLenVec.size(); i++) {
            qryLenVec.push_back(seqLenVec[i]);
        }
    }

    return make_tuple(seqLenVec[0], qryLenVec);
}


/**
 * Save the result of failCall in recall
 *
 * @param chrStartLenInfoGtTupMap     Save the result of failCall in recall
 * @param outFileName                 Output file name
 * 
 * 
 * @return int             0
**/
int VCFRecall::saveFailCall(
    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap, 
    const string & outFileName
)
{
    SAVE SAVEClass(outFileName);

    stringstream outStream; // Use stringstream instead of string concatenation
    outStream.str().reserve(bufferSize_);

    // Iterate through the dictionary and save the vcf that is not found
    for (const auto& [_, StartLenInfoGtTupMap] : chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (const auto& [_, LenInfoGtTup] : StartLenInfoGtTupMap) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            outStream << get<2>(LenInfoGtTup) + "\n";

            // save result
            if (outStream.tellp() >= bufferSize_) {  // Cache size is 10mb
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                // Clear stringstream
                outStream.str(string());
                outStream.clear();
            }
        }
    }

    // last save
    if (outStream.tellp() > 0) {  // Cache size is 10mb
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        // Clear stringstream
        outStream.str(string());
        outStream.clear();
    }

    return 0;
}


/**
 * Save the result of failCall in recall
 *
 * @param allLengthList     All variation lengths
 * 
 * 
 * @return int             0
**/
void VCFRecall::roc_calculate(
    const vector<int64_t> & allLengthList
) {
    // save result
    ofstream allFile;
    allFile.open("weight.all.table", ios::out);
    ofstream snpFile;
    snpFile.open("weight.snp.table", ios::out);
    ofstream indelFile;
    indelFile.open("weight.indel.table", ios::out);
    ofstream delFile;
    delFile.open("weight.del.table", ios::out);
    ofstream insFile;
    insFile.open("weight.ins.table", ios::out);
    
    // Calculate all the lengths to calculate roc
    // Cycle from big to small
    string outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    int64_t allTruePositives = allLengthList.size(); // All correct number of sites
    int64_t truePositives = 0; // True positive
    int64_t falsePositives = 0; // False positive
    int64_t trueNegatives = 0; // True negative
    int64_t falseNegatives = 0; // false negative

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 == rocRecallMap_.end()) {
            falsePositives += iter1->second.size();
        }
        else {
            truePositives += findIter1->second.size();
            falsePositives += iter1->second.size() - findIter1->second.size();
        }
        falseNegatives = allTruePositives-truePositives;
        float recall = (float)truePositives / allTruePositives;
        float precision = (float)truePositives / (truePositives + falsePositives);
        float Fscore = (2*recall*precision)/(recall + precision);

        outRoc += to_string(iter1->first) + "\t" 
                + to_string(truePositives) + "\t" 
                + to_string(falsePositives) + "\t" 
                + to_string(trueNegatives) + "\t" 
                + to_string(falseNegatives) + "\t"
                + to_string(recall) + "\t"
                + to_string(precision) + "\t"
                + to_string(Fscore) + "\n";
    }
    allFile << outRoc;
    allFile.close();


    // snp
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0; // All correct number of sites
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] == 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // True positive
    falsePositives = 0; // False positive
    trueNegatives = 0; // True negative
    falseNegatives = 0; // false negative

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint32_t allPositiveTmp = 0;
        // Number of SNPS in all
        // Traverse the length vector under the score
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] == 0) {
                allPositiveTmp++;
            }
        }

        uint32_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // Traverse the length vector under the score
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] == 0) {
                    truePositiveTmp++;
                }
            }
        }
        // If the number is 0, skip the score
        if (allPositiveTmp == 0) {
            continue;
        }

        truePositives += truePositiveTmp;
        falsePositives += allPositiveTmp - truePositiveTmp;

        falseNegatives = allTruePositives-truePositives;
        float recall = (float)truePositives / allTruePositives;
        float precision = (float)truePositives / (truePositives + falsePositives);
        float Fscore = (2*recall*precision)/(recall + precision);

        outRoc += to_string(iter1->first) + "\t" 
                + to_string(truePositives) + "\t" 
                + to_string(falsePositives) + "\t" 
                + to_string(trueNegatives) + "\t" 
                + to_string(falseNegatives) + "\t"
                + to_string(recall) + "\t"
                + to_string(precision) + "\t"
                + to_string(Fscore) + "\n";
    }
    snpFile << outRoc;
    snpFile.close();


    // indel
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0; // All correct number of sites
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] > -50 && allLengthList[i] < 50 && allLengthList[i] != 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // True positive
    falsePositives = 0; // False positive
    trueNegatives = 0; // True negative
    falseNegatives = 0; // false negative

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // Number of SNPS in all
        // Traverse the length vector under the score
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] > -50 && iter1->second[i] < 50 && iter1->second[i] != 0) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // Traverse the length vector under the score
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0) {
                    truePositiveTmp++;
                }
            }
        }
        // If the number is 0, skip the score
        if (allPositiveTmp == 0) {
            continue;
        }

        truePositives += truePositiveTmp;
        falsePositives += allPositiveTmp - truePositiveTmp;

        falseNegatives = allTruePositives-truePositives;
        float recall = (float)truePositives / allTruePositives;
        float precision = (float)truePositives / (truePositives + falsePositives);
        float Fscore = (2*recall*precision)/(recall + precision);

        outRoc += to_string(iter1->first) + "\t" 
                + to_string(truePositives) + "\t" 
                + to_string(falsePositives) + "\t" 
                + to_string(trueNegatives) + "\t" 
                + to_string(falseNegatives) + "\t"
                + to_string(recall) + "\t"
                + to_string(precision) + "\t"
                + to_string(Fscore) + "\n";
    }
    indelFile << outRoc;
    indelFile.close();


    // del
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0; // All correct number of sites
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] <= -50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // True positive
    falsePositives = 0; // False positive
    trueNegatives = 0; // True negative
    falseNegatives = 0; // false negative

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // Number of SNPS in all
        // Traverse the length vector under the score
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] <= -50) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // Traverse the length vector under the score
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] <= -50) {
                    truePositiveTmp++;
                }
            }
        }
        // If the number is 0, skip the score
        if (allPositiveTmp == 0) {
            continue;
        }

        truePositives += truePositiveTmp;
        falsePositives += allPositiveTmp - truePositiveTmp;

        falseNegatives = allTruePositives-truePositives;
        float recall = (float)truePositives / allTruePositives;
        float precision = (float)truePositives / (truePositives + falsePositives);
        float Fscore = (2*recall*precision)/(recall + precision);

        outRoc += to_string(iter1->first) + "\t" 
                + to_string(truePositives) + "\t" 
                + to_string(falsePositives) + "\t" 
                + to_string(trueNegatives) + "\t" 
                + to_string(falseNegatives) + "\t"
                + to_string(recall) + "\t"
                + to_string(precision) + "\t"
                + to_string(Fscore) + "\n";
    }
    delFile << outRoc;
    delFile.close();


    // ins
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0; // All correct number of sites
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] >= 50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // True positive
    falsePositives = 0; // False positive
    trueNegatives = 0; // True negative
    falseNegatives = 0; // false negative

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int64_t allPositiveTmp = 0;
        // Number of SNPS in all
        // Traverse the length vector under the score
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] >= 50) {
                allPositiveTmp++;
            }
        }

        int64_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // Traverse the length vector under the score
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] >= 50) {
                    truePositiveTmp++;
                }
            }
        }
        // If the number is 0, skip the score
        if (allPositiveTmp == 0) {
            continue;
        }

        truePositives += truePositiveTmp;
        falsePositives += allPositiveTmp - truePositiveTmp;

        falseNegatives = allTruePositives-truePositives;
        float recall = (float)truePositives / allTruePositives;
        float precision = (float)truePositives / (truePositives + falsePositives);
        float Fscore = (2*recall*precision)/(recall + precision);

        outRoc += to_string(iter1->first) + "\t" 
                + to_string(truePositives) + "\t" 
                + to_string(falsePositives) + "\t" 
                + to_string(trueNegatives) + "\t" 
                + to_string(falseNegatives) + "\t"
                + to_string(recall) + "\t"
                + to_string(precision) + "\t"
                + to_string(Fscore) + "\n";
    }
    insFile << outRoc;
    insFile.close();
}