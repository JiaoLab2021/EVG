#include "../include/vcf_recall.hpp"

#include <exception>
#include <regex>
#include <stdexcept>

#include "../include/logger.hpp"
#include "../include/vcf_sort.hpp"

using namespace std;

// Store the inputs, the thread count, the length bins and the ROC field (column 7 for INFO, 8 for FORMAT; INFO keys get their regex here).
VCFRecall::VCFRecall(
    const string& trueVCFFileName,
    const string& evaluateFileName,
    const string& rocKey,
    const string& rocType,
    uint32_t threads
) : trueVCFFileName_(trueVCFFileName), evaluateFileName_(evaluateFileName), rocKey_(rocKey), rocType_(rocType), threads_(threads)
{
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
    lengthVec_.push_back("-49/-20");
    lengthVec_.push_back("-19/-1");
    lengthVec_.push_back("0/0");
    lengthVec_.push_back("1/19");
    lengthVec_.push_back("20/49");
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

    rocColNum_ = 0;
    if (!rocKey_.empty())
    {
        rocKeyVec_ = split(rocKey_, ".");

        if (rocKeyVec_[0] == "INFO")
        {
            rocColNum_ = 7;
        }
        else if (rocKeyVec_[0] == "FORMAT")
        {
            rocColNum_ = 8;
        }

        if (rocKeyVec_.size() >= 2 && rocColNum_ != 8) {
            rocInfoReg_ = regex(rocKeyVec_[1] + "=(\\d+)");
        }
    }
}

// Index the non-reference records of the truth VCF by chromosome and position, checking its order in the same pass.
void VCFRecall::build_true_index()
{
    log_stream() << "Building index ...\n";

    SortedIndexPass pass(trueVCFFileName_);
    try {
        VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN VCFOPENCLASS(trueVCFFileName_, threads_);

        while (VCFOPENCLASS.read(INFOSTRUCTTMP)) {
            if (INFOSTRUCTTMP.line.empty() || INFOSTRUCTTMP.line.find("#") != string::npos) {
                continue;
            }

            pass.check_order(INFOSTRUCTTMP);
            pass.index([&] { index_true_record(INFOSTRUCTTMP); });
        }
    } catch (...) {
        pass.fail<bool>(current_exception()).get();
    }

    pass.finish(true).get();
}

// Add one truth record (skipped with a warning when its genotype has no ALT allele) to the truth index.
void VCFRecall::index_true_record(const VCFINFOSTRUCT & INFOSTRUCTTMP)
{
    vector<int> gtVec = get_gt(
        INFOSTRUCTTMP.lineVec
    );
    if (gtVec.empty() || (all_of(gtVec.begin(), gtVec.end(), [](int c) { return c == 0; }))) {
        warning_stream() << "The genotype of the variant is all '0' or all '.', skip this site -> "
            << INFOSTRUCTTMP.CHROM << " "
            << INFOSTRUCTTMP.POS << endl;
        return;
    }

    uint32_t refLen;
    vector<uint32_t> qryLenVec;
    tie(refLen, qryLenVec) = get_hap_len(
        INFOSTRUCTTMP.ID,
        INFOSTRUCTTMP.REF,
        INFOSTRUCTTMP.ALT,
        gtVec,
        "hap"
    );

    int32_t svLength = sv_length_select(
        refLen,
        qryLenVec,
        gtVec
    );

    trueVCFStrcuture_.allLengthList.push_back(svLength);

    trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);

    trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][INFOSTRUCTTMP.POS] = make_tuple(refLen, qryLenVec, INFOSTRUCTTMP.line, svLength, gtVec);
}

// Rewrite gramtools output with diploid genotypes from COV into convert.<name> and evaluate that file instead.
void VCFRecall::convert_gramtools()
{
    string outInformations;

    check_vcf_sort(evaluateFileName_);

    vector<string> prefixVec = split(evaluateFileName_, "/");
    string outFileName = "convert." + prefixVec.back();
    SAVE SaveClass(outFileName, threads_);

    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_, threads_);
    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            SaveClass.save(INFOSTRUCTTMP.line);
        }

        vector<string> gtVector = split(INFOSTRUCTTMP.lineVec.back(), ":");
        vector<string> formatVector = split(INFOSTRUCTTMP.FORMAT, ":");
        string gt;
        if (gtVector.size() > 0) {
            gt = gtVector[0];
        } else {
            gt = "0";
        }

        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
        int covIndex = 0;
        if (covItera != formatVector.end()) {
            covIndex = distance(formatVector.begin(), covItera);
        } else {
            outInformations += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        for (size_t i = 0; i < INFOSTRUCTTMP.lineVec.size()-1; i++) {
            if (i == 6) {
                outInformations += "PASS\t";
            } else {
                outInformations += INFOSTRUCTTMP.lineVec[i] + "\t";
            }
        }

        string siteCoverage = gtVector[covIndex];

        if (gt == "0") {
            outInformations += "\t0/0";
        } else {
            if (siteCoverage.find("0,") < siteCoverage.length()) {
                outInformations += "\t1/1";
            } else {
                outInformations += "\t0/1";
            }
        }

        for (size_t i = 1; i < gtVector.size(); i++) {
            outInformations += ":" + gtVector[i];
        }

        outInformations += "\n";

        if (outInformations.size() > bufferSize_) {
            SaveClass.save(outInformations);

            outInformations.clear();
        }
    }

    SaveClass.save(outInformations);

    evaluateFileName_ = outFileName;
}

// Evaluate calls by exact position: write the true, wrong-genotype, miscall and failcall VCFs, vcf_evulate.out and the ROC tables.
void VCFRecall::evaluate_genotype()
{
    log_stream() << "Evaluating the genotype results ..." << endl;

    check_vcf_sort(evaluateFileName_);

    // ------ Output files ------
    stringstream truetxtStream;
    truetxtStream.str().reserve(bufferSize_);
    SAVE trueFile("genotype.true.vcf.gz", threads_);

    stringstream genotypeMisTxtStream;
    genotypeMisTxtStream.str().reserve(bufferSize_);
    SAVE genotypeMisFile("genotype.err.vcf.gz", threads_);

    stringstream misCallTxtStream;
    misCallTxtStream.str().reserve(bufferSize_);
    SAVE misCallFile("miscall.vcf.gz", threads_);

    stringstream failCallTxtStream;
    failCallTxtStream.str().reserve(bufferSize_);
    SAVE failCallFile("failcall.vcf.gz", threads_);

    vector<int64_t> genotypeLenVec;
    vector<int64_t> misgenotypeLenVec;
    vector<int64_t> miscallLenVec;

    // ------ Classify each PASS call: true genotype, wrong genotype or miscall ------
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_, threads_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        if (INFOSTRUCTTMP.FILTER != "PASS") {
            warning_stream() << "'" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << " "
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        if (gtVec.empty() || (all_of(gtVec.begin(), gtVec.end(), [](int c) { return c == 0; }))) {
            warning_stream() << "The genotype of the variant is all '0' or all '.', skip this site -> "
                << INFOSTRUCTTMP.CHROM << " "
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        float rocNum = 0.0;
        rocNum = get_roc(INFOSTRUCTTMP);

        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID,
            INFOSTRUCTTMP.REF,
            INFOSTRUCTTMP.ALT,
            gtVec,
            "hap"
        );

        int32_t svLength = sv_length_select(
            refLen,
            qryLenVec,
            gtVec
        );

        rocCallMap_[rocNum].push_back(svLength);

        // ------ Match the alleles of the truth record at the same position ------
        auto iter1 = trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM);
        if (iter1 == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            warning_stream() << "" << INFOSTRUCTTMP.CHROM << " not in true set.\n";
            miscallLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }

        uint32_t genotypeTrueNum = 0;
        uint32_t trueRefLen;
        string trueVcfInfo;
        int32_t trueSvLen;

        auto iter2 = iter1->second.find(INFOSTRUCTTMP.POS);
        if (iter2 != iter1->second.end()) {
            trueRefLen = get<0>(iter2->second);
            const vector<uint32_t>& trueQryLenVec = get<1>(iter2->second);
            trueVcfInfo = get<2>(iter2->second);
            trueSvLen = get<3>(iter2->second);
            vector<int> trueGtVec = get<4>(iter2->second);

            if (trueQryLenVec.size() == 0) {
                throw runtime_error("trueQryLenVec.size() == 0 -> " + INFOSTRUCTTMP.CHROM + " " + to_string(INFOSTRUCTTMP.POS));
            }

            int64_t trueQryLenVecIdx = -1;
            for (size_t i = 0; i < qryLenVec.size(); i++) {
                uint32_t qryLen = qryLenVec[i];
                int callGt = gtVec[i];

                for (size_t j = 0; j < trueQryLenVec.size(); j++) {
                    if (static_cast<int64_t>(j) == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[j];
                    int trueGt = trueGtVec[j];

                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    if (refLen >= 50 && qryLen < 50) {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {
                        if ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {
                        if (((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {
                        if (((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen <= 2 && qryLen >= 3) {
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
                stop:;
            }
        }

        // ------ Record the class and flush every 10 MB ------
        if (genotypeTrueNum == 0) {
            miscallLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        } else {
            if (genotypeTrueNum >= gtVec.size()) {
                genotypeLenVec.push_back(trueSvLen);
                truetxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" +
                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                rocRecallMap_[rocNum].push_back(svLength);
            } else {
                misgenotypeLenVec.push_back(trueSvLen);
                genotypeMisTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" +
                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                if (rocType_ == "recall") {
                    rocRecallMap_[rocNum].push_back(svLength);
                }
            }
            trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
        }

        if (truetxtStream.tellp() >= bufferSize_ ||
            genotypeMisTxtStream.tellp() >= bufferSize_ ||
            misCallTxtStream.tellp() >= bufferSize_)
        {
            string truetxt = truetxtStream.str();
            trueFile.save(truetxt);
            truetxtStream.str(string());
            truetxtStream.clear();

            string genotypeMisTxt = genotypeMisTxtStream.str();
            genotypeMisFile.save(genotypeMisTxt);
            genotypeMisTxtStream.str(string());
            genotypeMisTxtStream.clear();

            string misCallTxt = misCallTxtStream.str();
            misCallFile.save(misCallTxt);
            misCallTxtStream.str(string());
            misCallTxtStream.clear();
        }
    }

    if (truetxtStream.tellp() > 0 ||
        genotypeMisTxtStream.tellp() > 0 ||
        misCallTxtStream.tellp() > 0)
    {
        string truetxt = truetxtStream.str();
        trueFile.save(truetxt);
        truetxtStream.str(string());
        truetxtStream.clear();

        string genotypeMisTxt = genotypeMisTxtStream.str();
        genotypeMisFile.save(genotypeMisTxt);
        genotypeMisTxtStream.str(string());
        genotypeMisTxtStream.clear();

        string misCallTxt = misCallTxtStream.str();
        misCallFile.save(misCallTxt);
        misCallTxtStream.str(string());
        misCallTxtStream.clear();
    }

    // ------ Truth records never matched ------
    for (auto it1 : trueVCFStrcuture_.chrStartLenInfoGtTupMap) {
        for (auto it2 : it1.second) {
            failCallTxtStream << get<2>(it2.second) + "\n";

            if (failCallTxtStream.tellp() >= bufferSize_) {
                string failCallTxt = failCallTxtStream.str();
                failCallFile.save(failCallTxt);
                failCallTxtStream.str(string());
                failCallTxtStream.clear();
            }
        }
    }
    if (failCallTxtStream.tellp() > 0) {
        string failCallTxt = failCallTxtStream.str();
        failCallFile.save(failCallTxt);
        failCallTxtStream.str(string());
        failCallTxtStream.clear();
    }

    // ------ Summary in vcf_evulate.out (appended), then ROC tables ------
    int64_t sv_genotype_recall = genotypeLenVec.size();
    int64_t sv_misgenotype_recall = misgenotypeLenVec.size();
    int64_t sv_mis_call = miscallLenVec.size();
    int64_t sv_all = trueVCFStrcuture_.allLengthList.size();

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

    for (size_t i = 0; i < lengthVec_.size(); i++) {
        outFile << lengthVec_[i] << ": "
                << genotype_length_count[i] << "/"
                << misgenotype_length_count[i] << "/"
                << genotype_length_count[i] + misgenotype_length_count[i] << "/"
                << mis_call_length_count[i] << "/"
                << genotype_length_count[i] + misgenotype_length_count[i] + mis_call_length_count[i] << "/"
                << all_length_count[i] - genotype_length_count[i] - misgenotype_length_count[i] << "/"
                << all_length_count[i] << endl;

        if (10 <= i && i <= 14) {
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

    outFile.close();

    roc_calculate(trueVCFStrcuture_.allLengthList);
}

// Evaluate calls by nearby position and length tolerance: write the recall VCFs, vcf_evulate.out and the ROC tables.
void VCFRecall::evaluate_recall()
{
    log_stream() << "Evaluating the genotype results ..." << endl;

    check_vcf_sort(evaluateFileName_);

    // ------ Output files ------
    stringstream truetxtStream;
    truetxtStream.str().reserve(bufferSize_);
    SAVE trueFile("recall.true.vcf.gz", threads_);

    stringstream trueGtTxtStream;
    trueGtTxtStream.str().reserve(bufferSize_);
    SAVE trueGtFile("genotype.true.vcf.gz", threads_);

    stringstream genotypeMisTxtStream;
    genotypeMisTxtStream.str().reserve(bufferSize_);
    SAVE genotypeMisFile("genotype.err.vcf.gz", threads_);

    stringstream misCallTxtStream;
    misCallTxtStream.str().reserve(bufferSize_);
    SAVE misCallFile("miscall.vcf.gz", threads_);

    vector<int64_t> genotypeLenVec;
    vector<int64_t> callLenVec;
    vector<int64_t> recallLenVec;

    string preChr;

    int leftIdxTmp = 0;
    int rightIdxTmp = 0;

    // ------ Classify each PASS call: recalled (true or wrong genotype) or miscall ------
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_, threads_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        if (INFOSTRUCTTMP.FILTER != "PASS") {
            warning_stream() << "'" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << " "
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        if (gtVec.empty() || (all_of(gtVec.begin(), gtVec.end(), [](int c) { return c == 0; }))) {
            warning_stream() << "The genotype of the variant is all '0' or all '.', skip this site -> "
                << INFOSTRUCTTMP.CHROM << " "
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        float rocNum = 0.0;
        rocNum = get_roc(INFOSTRUCTTMP);

        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID,
            INFOSTRUCTTMP.REF,
            INFOSTRUCTTMP.ALT,
            gtVec,
            "hap"
        );

        int32_t svLength = sv_length_select(
            refLen,
            qryLenVec,
            gtVec
        );

        rocCallMap_[rocNum].push_back(svLength);

        if (trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM) == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            warning_stream() << "'" << INFOSTRUCTTMP.CHROM << "' is not in the true set.\n";
            callLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // ------ Truth records within 10 bp (both alleles < 50 bp) or 200 bp, matched by class and length ------
        if (preChr != INFOSTRUCTTMP.CHROM) {
            leftIdxTmp = 0;
            rightIdxTmp = 0;
        }

        int64_t indexLeft = -1;
        int64_t indexRight = -1;

        uint32_t genotypeTrueNum = 0;

        uint32_t trueRefStart;
        uint32_t trueRefLen;
        int32_t trueSvLen = 0;
        string truevcfInfo;

        int64_t trueQryLenVecIdx = -1;

        for (size_t i = 0; i < qryLenVec.size(); i++) {
            uint32_t qryLen = qryLenVec[i];
            int callGt = gtVec[i];

            if (indexLeft == -1 && indexRight == -1) {
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
                    throw runtime_error("out of index, please check the data or code.");
                }
            }

            for (int64_t j = indexLeft; j <= indexRight; j++) {
                trueRefStart = trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM][j];

                if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(trueRefStart) == trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                    continue;
                }

                trueRefLen = get<0>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                const vector<uint32_t>& trueQryLenVec = get<1>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                truevcfInfo = get<2>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                trueSvLen = get<3>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                vector<int> trueGtVec = get<4>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);

                for (size_t k = 0; k < trueQryLenVec.size(); k++) {
                    if (static_cast<int64_t>(k) == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[k];
                    int trueGt = trueGtVec[k];

                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    if (refLen >= 50 && qryLen < 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }  else if (refLen <= 2 && qryLen >= 3) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            genotypeTrueNum++;

                            indexLeft = j;
                            indexRight = j;

                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }
                }
            }
            stop:;
        }

        // ------ Record the class and flush every 10 MB ------
        if (genotypeTrueNum > 0) {
            callLenVec.push_back(trueSvLen);

            truetxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" +
                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
            recallLenVec.push_back(trueSvLen);

            if (genotypeTrueNum >= gtVec.size()) {
                genotypeLenVec.push_back(trueSvLen);
                trueGtTxtStream << "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" +
                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                rocRecallMap_[rocNum].push_back(svLength);
            } else {
                genotypeMisTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" +
                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                if (rocType_ == "recall")
                {
                    rocRecallMap_[rocNum].push_back(svLength);
                }
            }

            if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(INFOSTRUCTTMP.POS) != trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
            }
        } else {
            callLenVec.push_back(svLength);
            misCallTxtStream << "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
        }

        if (truetxtStream.tellp() >= bufferSize_ ||
            trueGtTxtStream.tellp() >= bufferSize_ ||
            genotypeMisTxtStream.tellp() >= bufferSize_ ||
            misCallTxtStream.tellp() >= bufferSize_)
        {
            string truetxt = truetxtStream.str();
            trueFile.save(truetxt);
            truetxtStream.str(string());
            truetxtStream.clear();

            string trueGtTxt = trueGtTxtStream.str();
            trueGtFile.save(trueGtTxt);
            trueGtTxtStream.str(string());
            trueGtTxtStream.clear();

            string genotypeMisTxt = genotypeMisTxtStream.str();
            genotypeMisFile.save(genotypeMisTxt);
            genotypeMisTxtStream.str(string());
            genotypeMisTxtStream.clear();

            string misCallTxt = misCallTxtStream.str();
            misCallFile.save(misCallTxt);
            misCallTxtStream.str(string());
            misCallTxtStream.clear();
        }
    }

    if (truetxtStream.tellp() > 0 ||
        trueGtTxtStream.tellp() > 0 ||
        genotypeMisTxtStream.tellp() > 0 ||
        misCallTxtStream.tellp() > 0)
    {
        string truetxt = truetxtStream.str();
        trueFile.save(truetxt);
        truetxtStream.str(string());
        truetxtStream.clear();

        string trueGtTxt = trueGtTxtStream.str();
        trueGtFile.save(trueGtTxt);
        trueGtTxtStream.str(string());
        trueGtTxtStream.clear();

        string genotypeMisTxt = genotypeMisTxtStream.str();
        genotypeMisFile.save(genotypeMisTxt);
        genotypeMisTxtStream.str(string());
        genotypeMisTxtStream.clear();

        string misCallTxt = misCallTxtStream.str();
        misCallFile.save(misCallTxt);
        misCallTxtStream.str(string());
        misCallTxtStream.clear();
    }

    // ------ Truth records never matched ------
    saveFailCall(
        trueVCFStrcuture_.chrStartLenInfoGtTupMap,
        "failcall.vcf.gz"
    );

    // ------ Summary in vcf_evulate.out (appended), then ROC tables ------
    int64_t all_num = trueVCFStrcuture_.allLengthList.size();
    int64_t call_num = callLenVec.size();
    int64_t recall_num = recallLenVec.size();

    int64_t genotype_recall_number = genotypeLenVec.size();

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

    for (size_t i = 0; i < lengthVec_.size(); i++) {
        outFile << lengthVec_[i] << ": "
                << genotypeRecallNumVec[i] << "/"
                << recallNumVec[i] - genotypeRecallNumVec[i] << "/"
                << recallNumVec[i] << "/"
                << callNumVec[i] - recallNumVec[i] << "/"
                << callNumVec[i] << "/"
                << allNumVec[i] - recallNumVec[i] << "/"
                << allNumVec[i] << endl;

        if (10 <= i && i <= 14) {
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

    outFile.close();

    roc_calculate(trueVCFStrcuture_.allLengthList);
}

// Return the genotype of one sample (0 = last column) as integers; '.' counts as 0, unparsable genotypes become {0, 0}.
vector<int> VCFRecall::get_gt(
    const vector<string> & lineVec,
    int sampleIdx
) const {
    vector <int> gtVec;

    vector<string> formatVec = split(lineVec[8], ":");

    uint32_t gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));

    if (gtIndex == formatVec.size()) {
        warning_stream() << "No [GT] information in FORMAT, replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {
        string gt;

        if (sampleIdx == 0) {
            gt = VCFOPEN::sample_gt(lineVec, lineVec.size() - 1, gtIndex);
        } else {
            gt = VCFOPEN::sample_gt(lineVec, sampleIdx, gtIndex);
        }

        string splitStr;
        if (gt.find("/") != string::npos) {
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {
            splitStr = "|";
        } else {
            warning_stream() << "Genotype does not contain '\\' or '|', replace with 0|0 -> " << lineVec[0] << ":" << lineVec[1] << endl;
            gtVec = {0, 0};
            return gtVec;
        }

        auto splitResult = split(gt, splitStr);
        for (auto it : splitResult) {
            if (it == ".") {
                it = "0";
            }
            try {
                gtVec.push_back(stoi(it));
            } catch (const std::exception& e) {
                warning_stream() << "stof(it) conversion failed. Replaced with {0,0} -> " << it << endl;
                gtVec = {0, 0};
                return gtVec;
            }
        }
    }

    return gtVec;
}

// Return the ROC score of a record from the -s field (0 without -s); throw if the field is missing.
float VCFRecall::get_roc(
    const VCFINFOSTRUCT& INFOSTRUCTTMP
)
{
    if (rocKeyVec_.size() < 2) {
        return 0.0;
    }

    if (rocColNum_ == 8) {
        vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec.back(), ":");

        vector<string> rocInfVec = split(INFOSTRUCTTMP.FORMAT, ":");
        vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec_[1]);

        if (rocItera == rocInfVec.end()) {
            throw runtime_error("" + rocKeyVec_[1] + " not in " + rocKeyVec_[0] + " column -> " + INFOSTRUCTTMP.line);
        }

        return stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
    } else {
        smatch patternResult;
        if (!regex_search(INFOSTRUCTTMP.INFO, patternResult, rocInfoReg_)) {
            throw runtime_error("" + rocKeyVec_[1] + " not in " + rocKeyVec_[0] + " column -> " + INFOSTRUCTTMP.line);
        }

        string rocNumString = patternResult[1];

        return stof(rocNumString);
    }
}

// Return the length difference between the last non-reference allele of gtVec and REF.
int32_t VCFRecall::sv_length_select(
    const uint32_t & refLen,
    const vector<uint32_t> & qryLenVec,
    const vector<int> & gtVec
) const {
    int32_t svLength = 0;

    for (size_t i = 0; i < gtVec.size(); i++) {
        if (gtVec[i] == 0) {
            continue;
        } else {
            svLength = qryLenVec[i] - refLen;
        }
    }

    return svLength;
}

// Count the lengths falling into each length bin.
vector<int64_t> VCFRecall::count_num(
    vector<int64_t> length_list
) {
    vector<int64_t> out_length;
    for (size_t i = 0; i < lengthVec_.size(); i++) {
        int64_t x1 = std::stoll(split(lengthVec_[i], "/")[0]);
        int64_t x2 = std::stoll(split(lengthVec_[i], "/")[1]);

        vector<int64_t>::size_type result = count_if(length_list.begin(), length_list.end(), f_mod_r(x1, x2));

        out_length.push_back(result);
    }
    return out_length;
}

// Write the truth records that no call matched.
int VCFRecall::saveFailCall(
    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap,
    const string & outFileName
)
{
    SAVE SAVEClass(outFileName, threads_);

    stringstream outStream;
    outStream.str().reserve(bufferSize_);

    for (const auto& [_, StartLenInfoGtTupMap] : chrStartLenInfoGtTupMap) {
        for (const auto& [_, LenInfoGtTup] : StartLenInfoGtTupMap) {
            outStream << get<2>(LenInfoGtTup) + "\n";

            if (outStream.tellp() >= bufferSize_) {
                string outTxt = outStream.str();
                SAVEClass.save(outTxt);
                outStream.str(string());
                outStream.clear();
            }
        }
    }

    if (outStream.tellp() > 0) {
        string outTxt = outStream.str();
        SAVEClass.save(outTxt);
        outStream.str(string());
        outStream.clear();
    }

    return 0;
}

// Write the cumulative ROC tables (all, SNP, indel, deletion, insertion) from the highest score down.
void VCFRecall::roc_calculate(
    const vector<int64_t> & allLengthList
) {
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

    // ------ All variants ------
    string outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    int64_t allTruePositives = allLengthList.size();
    int64_t truePositives = 0;
    int64_t falsePositives = 0;
    int64_t trueNegatives = 0;
    int64_t falseNegatives = 0;

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

    // ------ SNPs ------
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0;
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] == 0) {
            allTruePositives++;
        }
    }
    truePositives = 0;
    falsePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint32_t allPositiveTmp = 0;
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] == 0) {
                allPositiveTmp++;
            }
        }

        uint32_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] == 0) {
                    truePositiveTmp++;
                }
            }
        }
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

    // ------ Indels (1-49 bp) ------
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0;
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] > -50 && allLengthList[i] < 50 && allLengthList[i] != 0) {
            allTruePositives++;
        }
    }
    truePositives = 0;
    falsePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] > -50 && iter1->second[i] < 50 && iter1->second[i] != 0) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0) {
                    truePositiveTmp++;
                }
            }
        }
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

    // ------ Deletions (>= 50 bp) ------
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0;
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] <= -50) {
            allTruePositives++;
        }
    }
    truePositives = 0;
    falsePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] <= -50) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] <= -50) {
                    truePositiveTmp++;
                }
            }
        }
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

    // ------ Insertions (>= 50 bp) ------
    outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    allTruePositives = 0;
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] >= 50) {
            allTruePositives++;
        }
    }
    truePositives = 0;
    falsePositives = 0;
    trueNegatives = 0;
    falseNegatives = 0;

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int64_t allPositiveTmp = 0;
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] >= 50) {
                allPositiveTmp++;
            }
        }

        int64_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] >= 50) {
                    truePositiveTmp++;
                }
            }
        }
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
