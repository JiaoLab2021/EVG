// g++ vcf_recall.cpp -o vcf_recall -lz -O2

#include "../include/vcf_recall.hpp"

using namespace std;

int main_recall(int argc, char** argv)
{
    string trueVCFFileName;
    string evaluateFileName;

    string model = "recall";
    string rocKey = "";
    // roc�������
    string rocType = "genotype";

    // �������
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

    // �������Ƿ���ȷ
    if (trueVCFFileName.empty() || evaluateFileName.empty() || model.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_recall(argv);
        return 1;
    }

    // ���rocType�Ƿ���ȷ
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

// �����ĵ�
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

    bufferSize_ = 10 * 1024 * 1024;

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

        // ���ݻ����ͽ��й���
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //�����(0/0, .)��ʽ����ֱ�����������߷����˿��б�����������
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The genotype of the variant is empty, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // ��ȡref�͵����͵ĳ���
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // ��ȡ����ĳ���
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // �������б���ĳ���
        trueVCFStrcuture_.allLengthList.push_back(svLength);

        // ������Ϣ
        trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);

        // �ȳ�ʼ����ϣ��
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

    // check file�Ƿ�����
    check_vcf_sort(evaluateFileName_);

    // output file stream
    vector<string> prefixVec = split(evaluateFileName_, "/");
    string outFileName = "convert." + prefixVec[prefixVec.size()-1];
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

        // ��ȡinformationsVector�л����ͺ�λ�㸲�Ƕ���Ϣ
        vector<string> gtVector = split(INFOSTRUCTTMP.lineVec.back(), ":");
        vector<string> formatVector = split(INFOSTRUCTTMP.FORMAT, ":");
        string gt;
        if (gtVector.size() > 0) {
            gt = gtVector[0];
        } else {
            gt = "0";
        }

        // ��FORMAT�ֶ���COV��λ��
        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
        int covIndex = 0;
        if (covItera != formatVector.end()) { // FORMAT����COV
            covIndex = distance(formatVector.begin(), covItera);
        } else {  // û�еĻ�ֱ��ת����������һ��ѭ��
            outInformations += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // FILTER�ֶ�ǰ�ߵ����ݲ��䣬ֱ�Ӽӵ�outInformations��
        for (int i = 0; i < INFOSTRUCTTMP.lineVec.size()-1; i++) {
            if (i == 6) {  // FILTER�и�ΪPASS��Դ�ļ��е�Ϊ.
                outInformations += "PASS\t";
            } else {
                outInformations += INFOSTRUCTTMP.lineVec[i] + "\t";
            }
        }

        string siteCoverage = gtVector[covIndex];
        
        // ���gt=0�Ļ������λ��û�б��죬��gt��Ϊ0/0
        if (gt == "0") {
            outInformations += "\t0/0";
        } else {  // gt��Ϊ0��ʱ��
            // ���λ�㸲�Ƕ�����0,�ֶΣ������λ��Ϊ���͵�ͻ��λ��
            if (siteCoverage.find("0,") < siteCoverage.length()) {
                outInformations += "\t1/1";
            } else {  // ���λ�㸲�Ƕ�û��0,�ֶΣ������λ��Ϊ�Ӻϵ�ͻ��λ��
                outInformations += "\t0/1";
            }
        }

        // ��informations�����ֶμ���
        for (int i = 1; i < gtVector.size(); i++) {
            outInformations += ":" + gtVector[i];
        }

        outInformations += "\n";

        // ÿ10Mbд��һ���ļ�
        if (outInformations.size() > bufferSize_) {
            SaveClass.save(outInformations);

            // ����ַ���
            outInformations.clear();
        }
    }

    // ���д��һ���ļ�
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

    // check file�Ƿ�����
    check_vcf_sort(evaluateFileName_);
    
    // ���������ȷ�Ľ��
    string truetxt;
    SAVE trueFile("genotype.true.vcf.gz");

    // ������ʹ���Ľ��
    string genotypeMisTxt;
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // ����miscall�Ľ��
    string misCallTxt;
    SAVE misCallFile("miscall.vcf.gz");

    // û���ҵ��ı���
    string failCallTxt;
    SAVE failCallFile("failcall.vcf.gz");

    // ����vector
    vector<uint32_t> genotypeLenVec;
    vector<uint32_t> misgenotypeLenVec;
    vector<uint32_t> miscallLenVec;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // ����ע����
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }
        
        // ����FILTER�ֶν��й���
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // �������(PASS)��ֱ������
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // ��ȡ��������Ϣ
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // ���ݻ����ͽ��й���
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //�����(0/0, .)��ʽ����ֱ�����������߷����˿��б�����������
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

        // ��ȡref�͵����͵ĳ���
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // ��ȡ����ĳ���
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // ��ӳ��ȵ�rocAllMap��
        rocCallMap_[rocNum].push_back(svLength);

        // ���л����������
        // ���Ⱦɫ���Ƿ���ڣ������ھ���miscall
        auto iter1 = trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM);
        if (iter1 == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: " << INFOSTRUCTTMP.CHROM << " not in true set.\n";
            miscallLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // ��¼genotype�Ľ��
        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype
        uint32_t trueRefLen;
        string trueVcfInfo;
        int32_t trueSvLen;

        // �ڶ�ӦȾɫ���map�в�����ֵ��refLen qryLenVec
        auto iter2 = iter1->second.find(INFOSTRUCTTMP.POS);
        if (iter2 != iter1->second.end()) {  // ����ҵ��ˣ��������жϷ����Ƿ���ȷ
            // true������Ϣ
            trueRefLen = get<0>(iter2->second);
            const vector<uint32_t>& trueQryLenVec = get<1>(iter2->second);
            trueVcfInfo = get<2>(iter2->second);
            trueSvLen = get<3>(iter2->second);
            vector<int> trueGtVec = get<4>(iter2->second);

            // �������Ϊ0���򱨴����˳����롣
            if (trueQryLenVec.size() == 0) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: trueQryLenVec.size() == 0 -> " 
                    << INFOSTRUCTTMP.CHROM << " " 
                    << INFOSTRUCTTMP.POS << endl;
                exit(1);
            }

            // ������λ����
            int trueQryLenVecIdx = -1;  // ��¼��ֵ���ù��ĵ����ͣ���ֹ�ظ�����
            for (size_t i = 0; i < qryLenVec.size(); i++) {
                uint32_t qryLen = qryLenVec[i];
                int callGt = gtVec[i];

                for (size_t j = 0; j < trueQryLenVec.size(); j++) {
                    // �ù��ĵ���������
                    if (j == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[j];
                    int trueGt = trueGtVec[j];

                    // callGt��trueGt����һ����0������һ�����ǣ�����һ��ѭ������ֹSNP�ж�ʱ�����
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // ȱʧ
                    if (refLen >= 50 && qryLen < 50) {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // ����
                        if ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // �滻
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
                stop:; // ����ҵ��ˣ����˳�Ƕ��ѭ���������õ����͵�ѭ����������һ��������
            }
        }

        // �жϷ��͵Ľ��
        if (genotypeTrueNum == 0) {  // �ҵ��ı��첻���漯�� 
            miscallLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        } else {
            // ������ȷ
            if (genotypeTrueNum >= gtVec.size()) {
                genotypeLenVec.push_back(trueSvLen);
                truetxt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // ��ӳ��ȵ�rocTrueMap��
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // ���ʹ���
                misgenotypeLenVec.push_back(trueSvLen);
                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // ���ж�����recall������roc�������������ֻ��genotype����roc
                if (rocType_ == "recall") {
                    // ��ӳ��ȵ�rocTrueMap��
                    rocRecallMap_[rocNum].push_back(svLength);
                }
                
            }
            // ɾ���Ѿ����ù��ı���
            trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
        }

        // save result
        if (truetxt.size() > bufferSize_ || 
            genotypeMisTxt.size() > bufferSize_ || 
            misCallTxt.size() > bufferSize_) // ÿ10Mbд��һ��
        {
            trueFile.save(truetxt);  // ������ȷ
            genotypeMisFile.save(genotypeMisTxt);  // ���ʹ���
            misCallFile.save(misCallTxt);  // �����漯�еı���

            // ����ַ���
            truetxt.clear();
            genotypeMisTxt.clear();
            misCallTxt.clear();
        }
    }

    // save result
    trueFile.save(truetxt);  // ������ȷ
    genotypeMisFile.save(genotypeMisTxt);  // ���ʹ���
    misCallFile.save(misCallTxt);  // �����漯�еı���

    // δ�ҵ�����ʵ����
    for (auto it1 : trueVCFStrcuture_.chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (auto it2 : it1.second) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            failCallTxt += get<2>(it2.second) + "\n";

            if (failCallTxt.size() > bufferSize_) {  // ÿ10Mbд��һ��
                failCallFile.save(failCallTxt);  // δ�ҵ�����ʵ����

                // ����ַ���
                failCallTxt.clear();
            }
        }
    }
    failCallFile.save(failCallTxt);  // δ�ҵ�����ʵ����

    uint64_t sv_genotype_recall = genotypeLenVec.size();
    uint64_t sv_misgenotype_recall = misgenotypeLenVec.size();
    uint64_t sv_mis_call = miscallLenVec.size();
    uint64_t sv_all = trueVCFStrcuture_.allLengthList.size();

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

    vector<uint64_t> all_length_count = count_num(trueVCFStrcuture_.allLengthList);
    vector<uint64_t> mis_call_length_count = count_num(miscallLenVec);
    vector<uint64_t> genotype_length_count = count_num(genotypeLenVec);
    vector<uint64_t> misgenotype_length_count = count_num(misgenotypeLenVec);
    
    outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

    // ���Ƚ�����
    for (int i = 0; i < lengthVec_.size(); i++)
    {
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

    // �ͷ��ڴ�
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

    // ���vcf�ļ��Ƿ�����
    check_vcf_sort(evaluateFileName_);

    // ������ȷ�Ľ��
    string truetxt;
    SAVE trueFile("recall.true.vcf.gz");

    // ���������ȷ�Ľ��
    string true_Gt_txt;
    SAVE trueGtFile("genotype.true.vcf.gz");

    // ������ʹ���Ľ��
    string genotypeMisTxt;
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // ����miscall�Ľ��
    string misCallTxt;
    SAVE misCallFile("miscall.vcf.gz");

    // ����vector
    vector<uint32_t> genotypeLenVec;
    vector<uint32_t> callLenVec;
    vector<uint32_t> recallLenVec;

    string preChr; // �����ж��ǲ����µ�Ⱦɫ�壬�ǵĻ������¹���vector��

    // ��¼�Ѿ����ù��ı��죬��ֹ�ظ�����
    map<string,vector<uint32_t> > selectChrStartMap;

    // ���ֲ��ҷ���������
    int leftIdxTmp = 0;
    int rightIdxTmp = 0;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // ����ע����
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        // ����FILTER�ֶν��й���
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // �������(PASS)��ֱ������
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // ��ȡ��������Ϣ
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // ���ݻ����ͽ��й���
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //�����(0/0, .)��ʽ����ֱ�����������߷����˿��б�����������
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


        // ��ȡref�͵����͵ĳ���
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // ��ȡ����ĳ���
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // ��ӳ��ȵ�rocAllMap��
        rocCallMap_[rocNum].push_back(svLength);


        // recall
        // ���Ⱦɫ���Ƿ���ڣ������ھ���miscall
        if (trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM) == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.CHROM << "' is not in the true set.\n";
            callLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }
    

        if (preChr != INFOSTRUCTTMP.CHROM) {  // ���ֲ��ҵ������������㣬���ڼӿ��ѯ�ٶ�
            // ������������
            leftIdxTmp = 0;
            rightIdxTmp = 0;
        }

        // ��¼���ֲ��ҷ��������������λ������ͬһ��refStart
        int indexLeft = -1;
        int indexRight = -1;

        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype

        // ��¼recall��  trueRefStart, trueRefLen, trueSvLen��truevcfInfo
        uint32_t trueRefStart;
        uint32_t trueRefLen;
        int32_t trueSvLen;
        string truevcfInfo;

        // ��¼�漯�б��ù��ĵ�����
        int trueQryLenVecIdx = -1;

        // �������λ����
        for (size_t i = 0; i < qryLenVec.size(); i++) {
            uint32_t qryLen = qryLenVec[i];
            int callGt = gtVec[i];

            // ���ֲ��ҷ���200/10bp�ڵı���������
            if (indexLeft == -1 && indexRight == -1) {  // �µĵ�λ�����ٲ���
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
                } if (indexLeft < 0 || indexRight >= trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].size()) {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: out of index, please check the data or code.\n";
                    exit(1);
                }
            }

            // �������ֲ��ҷ�������
            for (int j = indexLeft; j <= indexRight; j++) {
                // true�������Ϣ
                trueRefStart = trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM][j];

                // �ȼ�������û�б��ù��������ɾ���˾���һ������
                if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(trueRefStart) == trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                    continue;
                }
                
                trueRefLen = get<0>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                const vector<uint32_t>& trueQryLenVec = get<1>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                truevcfInfo = get<2>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                trueSvLen = get<3>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                vector<int> trueGtVec = get<4>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);

                for (size_t k = 0; k < trueQryLenVec.size(); k++) {  // һ��λ���ж����λ�����ʱ��trueSeqLenVec�洢�˸�����λ����ĳ��ȣ���˱���������λ��ı����Ƿ��merge��ߵ�һ�£�����0
                    // �ù��ĵ���������
                    if (k == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[k];
                    int trueGt = trueGtVec[k];

                    // callGt��trueGt����һ����0������һ�����ǣ�����һ��ѭ������ֹSNP�ж�ʱ�����
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // ȱʧ
                    if (refLen >= 50 && qryLen < 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // ����
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // �滻
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {  // snp
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {  // indel-del
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }  else if (refLen <= 2 && qryLen >= 3) {  // indel-ins
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // ��¼������recall and genotype ��ȷ
                            genotypeTrueNum++;

                            // �����λ������ͬһ����ʼλ��
                            indexLeft = j;
                            indexRight = j;

                            // ��¼�ù������͵�����
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }
                }
            }
            stop:; // ����ҵ��ˣ����˳�Ƕ��ѭ���������õ����͵�ѭ����������һ��������
        }

        // �ж�Ѱ�ҵĽ������� 
        if (genotypeTrueNum > 0) {  // ����0�����ҵ���
            // ���call��vcf-vector
            callLenVec.push_back(trueSvLen);

            truetxt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
            recallLenVec.push_back(trueSvLen);

            // ������ȷ
            if (genotypeTrueNum >= gtVec.size()) {  // ���ڵ���gtVec.size()���������ȷ
                genotypeLenVec.push_back(trueSvLen);
                true_Gt_txt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                // ��ӳ��ȵ�rocTrueMap��
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // ���ʹ���
                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                // ���ж�����recall������roc�������������ֻ��genotype����roc
                if (rocType_ == "recall")
                {
                    // ��ӳ��ȵ�rocTrueMap��
                    rocRecallMap_[rocNum].push_back(svLength);
                }
            }

            // ɾ���Ѿ����ù��ı���
            if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(INFOSTRUCTTMP.POS) != trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
            }
        } else {  // ����ϱ�ѭ��û�ҵ��漯��ı��죬��������Լ��ҵ��ĳ��������
            callLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
        }

        // save result
        if (truetxt.size() > bufferSize_ || 
            true_Gt_txt.size() > bufferSize_ || 
            genotypeMisTxt.size() > bufferSize_ || 
            misCallTxt.size() > bufferSize_) // ÿ10Mbд��һ��
        {
            trueFile.save(truetxt);
            trueGtFile.save(true_Gt_txt);
            genotypeMisFile.save(genotypeMisTxt);  // ���ʹ���
            misCallFile.save(misCallTxt);  // �����漯�еı���

            // ����ַ���
            truetxt.clear();
            true_Gt_txt.clear();
            genotypeMisTxt.clear();
            misCallTxt.clear();
        }
    }

    // save result
    trueFile.save(truetxt);
    trueGtFile.save(true_Gt_txt);
    genotypeMisFile.save(genotypeMisTxt);  // ���ʹ���
    misCallFile.save(misCallTxt);  // �����漯�еı���

    // ����û���ҵ����죬������
    saveFailCall(
        trueVCFStrcuture_.chrStartLenInfoGtTupMap, 
        "failcall.vcf.gz"
    );

    uint64_t all_num = trueVCFStrcuture_.allLengthList.size();
    uint64_t call_num = callLenVec.size();
    uint64_t recall_num = recallLenVec.size();

    uint64_t genotype_recall_number = genotypeLenVec.size();

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

    vector<uint64_t> allNumVec = count_num(trueVCFStrcuture_.allLengthList);
    vector<uint64_t> callNumVec = count_num(callLenVec);
    vector<uint64_t> recallNumVec = count_num(recallLenVec);
    vector<uint64_t> genotypeRecallNumVec = count_num(genotypeLenVec);
    
    outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

    // ���Ƚ�����
    for (int i = 0; i < lengthVec_.size(); i++)
    {
        outFile << lengthVec_[i] << ": " 
                << genotypeRecallNumVec[i] << "/" 
                << recallNumVec[i] - genotypeRecallNumVec[i] << "/" 
                << recallNumVec[i] << "/" 
                << callNumVec[i] - recallNumVec[i] << "/"
                << callNumVec[i] << "/"
                << allNumVec[i] - recallNumVec[i] << "/"
                << allNumVec[i] << endl;
        
        if ( 10 <= i && i <= 12)
        {
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

    // �ͷ��ڴ�
    outFile.close();

    // ����ROC
    roc_calculate(trueVCFStrcuture_.allLengthList);
}


/**
 * ��ȡλ��������б�.
 *
 * @param lineVec          lineVec
 * @param sampleIdx        sample�����͵�����,Ĭ��ֵ0�������һ��
 * 
 * 
 * @return gtVec           vector <int>
**/
vector<int> VCFRecall::get_gt(
    const vector<string> & lineVec, 
    int sampleIdx
) {
    vector <int> gtVec;  // λ����͵�vector

    // FORMAT�ַ����
    vector<string> formatVec = split(lineVec[8], ":");

    int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // ��ȡGT������λ��

    if (gtIndex == formatVec.size()) {  // �ж�index�Ƿ���ڣ������ڵĻ����ػ����Ͷ�Ϊ0��
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {  // ������ڣ�����б���
        string gt;  // �洢�������ֶ�

        if (sampleIdx == 0) {  // û��ָ�������������һ��
            gt = split(lineVec[lineVec.size()-1], ":")[gtIndex];  // gt�ֶ�
        } else {
            gt = split(lineVec[sampleIdx], ":")[gtIndex];  // gt�ֶ�
        }
        
        string splitStr;  // gt�еķָ���
        if (gt.find("/") != string::npos) {  // �жϡ�/���ָ���
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {  // �жϡ�|��Ϊ�ָ���
            splitStr = "|";
        } else {  // ��֪����ʱ��Ϊ���ؿ�ֵ
            gtVec = {0, 0};
            return gtVec;
        }

        // �ҵ�gt�󣬶��䰴splitStr��ֲ�ѭ��
        auto splitResult = split(gt, splitStr);
        for (auto it : splitResult) {
            // ���Ϊ'.'��������λ��
            if (it == ".") {
                gtVec = {0, 0};
                return gtVec;
            }
            // Prevent error in stof(it) conversion
            try {
                gtVec.push_back(stoi(it));  // ��ӵ�vector��
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
 * ��ȡλ��������б�.
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

    if (rocColNum_ == 8) {  // FORMAT�ֶ�
        vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec.back(), ":");

        vector<string> rocInfVec = split(INFOSTRUCTTMP.FORMAT, ":");
        // rocNum���±�
        vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec_[1]);

        if (rocItera == rocInfVec.end()) {  // �����ֶ�����û�ж�Ӧ��rocNum��Ϣ
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: " << rocKeyVec_[1] << " not in " << rocKeyVec_[0] << " column -> " 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        return stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
    } else {  // INFO�ֶ�
        smatch patternResult; // ������ʽ�Ľ��
        regex pattern(rocKeyVec_[1] + "=(\\d+)");

        if (!regex_search(INFOSTRUCTTMP.INFO, patternResult, pattern)) {  // ����Ƿ�ƥ�䵽 rocKeyVec_[1] ��ֵ
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
 * ��ȡ����ĳ�����Ϣ
 *
 * @param refLen            ref����
 * @param qryLenVec         qry�����б�
 * @param gtVec             �������б�
 * 
 * 
 * @return int32_t          svLength
**/
int32_t VCFRecall::sv_length_select(
    const uint32_t & refLen, 
    const vector<uint32_t> & qryLenVec, 
    const vector<int> & gtVec
) {
    int32_t svLength;

    for (size_t i = 0; i < gtVec.size(); i++) {
        if (gtVec[i] == 0) {  // �����������0������
            continue;
        } else {
            svLength = qryLenVec[i] - refLen;
        }
    }
    
    return svLength;
}


/**
 * ͳ�Ʊ��쳤����Ϣ
 *
 * @param lengthVec_          ���ֵ�����
 * @param length_list        �����б�
 * 
 * 
 * @return vector<uint64_t>  ÿ������ĳ���
**/
vector<uint64_t> VCFRecall::count_num(
    vector<uint32_t> length_list
)
{
    vector<uint64_t> out_length;
    for (int i = 0; i < lengthVec_.size(); i++) {
        uint32_t x1 = std::stoul(split(lengthVec_[i], "/")[0]);
        uint32_t x2 = std::stoul(split(lengthVec_[i], "/")[1]);

        vector<uint64_t>::size_type result = count_if(length_list.begin(), length_list.end(), f_mod(x1, x2));

        out_length.push_back(result);
    }
    return out_length;
}


/**
 * ��ȡ�����Ͷ�Ӧ�ĳ�����Ϣ.
 *
 * @param svType                         ��������
 * @param refSeq                         ref����Ϣ
 * @param qrySeqs                        qry����Ϣ
 * @param gtVec                          λ��ķ�����Ϣ
 * @param lenType                        ֻȡ�����Ͷ�Ӧ�ĳ��Ȼ������еĳ���(hap/all)
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
    // ���ģʽ�Ƿ���ȷ������ȷ�˳�����
    if (lenType != "hap" && lenType != "all") {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "Error: lenType -> " 
            << lenType << endl;
        exit(1);
    }
    
    uint32_t refLen = refSeq.size();  // ref���г���
    vector<string> qrySeqVec = split(qrySeqs, ",");  // qry�������б�
    
    // ��������Ƿ�Խ�磬���ֻҪhap�ĳ���ʱ���ټ��
    if (lenType == "hap") {
        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
        if (maxGtNum > qrySeqVec.size()) {  // �ȼ���Ƿ������Ƿ�Խ��
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: number of genotyping and ALT sequences do not match -> " 
                << qrySeqs << endl;
            exit(1);
        }
    }

    // ����ref��qry��ͬ�ĳ�������
    vector<uint32_t> seqLenVec;
    seqLenVec.push_back(refLen);

    // ������λ�����б�
    for (size_t i = 0; i < qrySeqVec.size(); i++) {
        string qrySeq = qrySeqVec[i];

        // ��ʱ�洢����
        uint32_t refLenTmp = refLen;
        uint32_t qryLenTmp = qrySeq.size(); 

        if (qrySeq.find(">") != string::npos) {  // GraphTyper2�Ľ��
            // ʹ��������ʽ��ȡ������Ϣ
            std::regex reg("SVSIZE=(\\d+)");
            std::smatch match;
            // ���qry��������û�г�����Ϣ��û�еĻ�������λ��
            if (std::regex_search(qrySeq, match, reg)) {  // <INS:SVSIZE=97:BREAKPOINT1>
                // <INS:SVSIZE=90:BREAKPOINT1>
                if (qrySeq.find("<INS") != string::npos && qrySeq.find("<DEL") == string::npos) {  // �ַ����а�������
                    refLenTmp = refLen;
                    qryLenTmp = std::stoul(match[1].str());;
                } else if (qrySeq.find("<INS") == string::npos && qrySeq.find("<DEL") != string::npos) {  // �ַ����а���ȱʧ: <DEL:SVSIZE=1720:BREAKPOINT>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLen;
                } else if (qrySeq.find("<DUP") != string::npos) {  // �ַ����а����ظ����ֶΣ�GraphTyper2��<DUP:SVSIZE=2806:BREAKPOINT1>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLenTmp * 2;
                }
            } else {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No length information in column 4 -> " << qrySeq << endl;
                refLenTmp = refLen;
                qryLenTmp = qrySeq.size();
            }

            seqLenVec[0] = refLenTmp; // ����ref�ĳ���
        }

        // �ж�BayesTyper�Ľ��Ϊduplication����ref_lenΪ1-2��qry_seqȴ�ܳ�ʱ
        if (svType.find("Duplication") != string::npos && refLenTmp <= 2) {  // BayesTyper���duplication��ɲ��룬�������¼��㳤��
            refLenTmp = qryLenTmp;
            qryLenTmp *= 2;
        }

        // ���qry�ĳ���
        if (seqLenVec.size() == (i + 1)) {  // �жϵ�����ȷʵ���Լ���λ���ϣ�
            seqLenVec.push_back(qryLenTmp);  // ���qry�ĳ���
        } else {  // �������б��Ȳ�����ʱ����
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Warning: Incorrect index for haplotype length -> " << qrySeqs << endl;
            exit(1);
        }
    }

    // �һ����Ͷ�Ӧ��qry����
    vector<uint32_t> qryLenVec;
    if (lenType == "hap") {  // ֻȡ�����͵����г���
        if (gtVec.size() == 1) {  // gtֻ��һ��������Ϊ���ϵı��죬�������������һ���ĳ���
            if (gtVec[0] == 0) {  // ���������Ϊ0��������λ��
                // ���� '0/0'
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
    } else {  // ȡ���еĳ���
        for (size_t i = 1; i < seqLenVec.size(); i++) {
            qryLenVec.push_back(seqLenVec[i]);
        }
    }

    return make_tuple(seqLenVec[0], qryLenVec);
}


/**
 * ����recall��failCall�Ľ��
 *
 * @param chrStartLenInfoGtTupMap     �漯����ʣ�µ�vcf��Ϣ
 * @param outFileName                 ����ļ���
 * 
 * 
 * @return int             0
**/
int VCFRecall::saveFailCall(
    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap, 
    const string & outFileName
)
{
    // û���ҵ��ı���
    string failCallTxt;
    // ����ļ���
    gzFile gzfp = gzopen(outFileName.c_str(), "wb");

    // �����ֵ䣬��û���ҵ���vcf���б���
    for (const auto& [_, StartLenInfoGtTupMap] : chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (const auto& [_, LenInfoGtTup] : StartLenInfoGtTupMap) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            failCallTxt += get<2>(LenInfoGtTup) + "\n";

            // save result
            if (failCallTxt.size() > bufferSize_) {  // ÿ10Mbд��һ��
                gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

                // ����ַ���
                failCallTxt.clear();
            }
        }
    }
    // save result
    gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

    // �ͷ��ڴ棬�ر��ļ�
    gzclose(gzfp);

    return 0;
}


/**
 * ����recall��failCall�Ľ��
 *
 * @param allLengthList     ���б��쳤��
 * 
 * 
 * @return int             0
**/
void VCFRecall::roc_calculate(
    const vector<uint32_t> & allLengthList
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
    
    // �������еĳ��ȼ���roc
    // �Ӵ���Сѭ��
    string outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    uint64_t allTruePositives = allLengthList.size(); // ������ȷ��λ������
    uint64_t truePositives = 0; // ������
    uint64_t falsePositives = 0; // ������
    uint64_t trueNegatives = 0; // ������
    uint64_t falseNegatives = 0; // ������

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
    allTruePositives = 0; // ������ȷ��λ������
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] == 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // ������
    falsePositives = 0; // ������
    trueNegatives = 0; // ������
    falseNegatives = 0; // ������

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint32_t allPositiveTmp = 0;
        // all�е�snp����
        // ������score�µĳ���vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] == 0) {
                allPositiveTmp++;
            }
        }

        uint32_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // ������score�µĳ���vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] == 0) {
                    truePositiveTmp++;
                }
            }
        }
        // ���������0��������score
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
    allTruePositives = 0; // ������ȷ��λ������
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] > -50 && allLengthList[i] < 50 && allLengthList[i] != 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // ������
    falsePositives = 0; // ������
    trueNegatives = 0; // ������
    falseNegatives = 0; // ������

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // all�е�snp����
        // ������score�µĳ���vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] > -50 && iter1->second[i] < 50 && iter1->second[i] != 0) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // ������score�µĳ���vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0) {
                    truePositiveTmp++;
                }
            }
        }
        // ���������0��������score
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
    allTruePositives = 0; // ������ȷ��λ������
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] <= -50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // ������
    falsePositives = 0; // ������
    trueNegatives = 0; // ������
    falseNegatives = 0; // ������

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // all�е�snp����
        // ������score�µĳ���vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] <= -50) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // ������score�µĳ���vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] <= -50) {
                    truePositiveTmp++;
                }
            }
        }
        // ���������0��������score
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
    allTruePositives = 0; // ������ȷ��λ������
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] >= 50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // ������
    falsePositives = 0; // ������
    trueNegatives = 0; // ������
    falseNegatives = 0; // ������

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint64_t allPositiveTmp = 0;
        // all�е�snp����
        // ������score�µĳ���vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] >= 50) {
                allPositiveTmp++;
            }
        }

        uint64_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // ������score�µĳ���vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] >= 50) {
                    truePositiveTmp++;
                }
            }
        }
        // ���������0��������score
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