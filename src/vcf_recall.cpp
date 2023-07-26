// g++ vcf_recall.cpp -o vcf_recall -lz -O2

#include "../include/vcf_recall.hpp"

using namespace std;

int main_recall(int argc, char** argv)
{
    string trueVCFFileName;
    string evaluateFileName;

    string model = "recall";
    string rocKey = "";
    // roc计算规则
    string rocType = "genotype";

    // 输入参数
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

    // 检查参数是否正确
    if (trueVCFFileName.empty() || evaluateFileName.empty() || model.empty()) {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
             << "Parameter error.\n";
        help_recall(argv);
        return 1;
    }

    // 检查rocType是否正确
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

// 帮助文档
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

        // 根据基因型进行过滤
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: The genotype of the variant is empty, skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // 获取ref和单倍型的长度
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // 获取变异的长度
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // 保存所有变异的长度
        trueVCFStrcuture_.allLengthList.push_back(svLength);

        // 变异信息
        trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM].push_back(INFOSTRUCTTMP.POS);

        // 先初始化哈希表
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

    // check file是否排序
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

        // 提取informationsVector中基因型和位点覆盖度信息
        vector<string> gtVector = split(INFOSTRUCTTMP.lineVec.back(), ":");
        vector<string> formatVector = split(INFOSTRUCTTMP.FORMAT, ":");
        string gt;
        if (gtVector.size() > 0) {
            gt = gtVector[0];
        } else {
            gt = "0";
        }

        // 找FORMAT字段中COV的位置
        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
        int covIndex = 0;
        if (covItera != formatVector.end()) { // FORMAT中有COV
            covIndex = distance(formatVector.begin(), covItera);
        } else {  // 没有的话直接转换，继续下一个循环
            outInformations += INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // FILTER字段前边的内容不变，直接加到outInformations上
        for (int i = 0; i < INFOSTRUCTTMP.lineVec.size()-1; i++) {
            if (i == 6) {  // FILTER列改为PASS，源文件有的为.
                outInformations += "PASS\t";
            } else {
                outInformations += INFOSTRUCTTMP.lineVec[i] + "\t";
            }
        }

        string siteCoverage = gtVector[covIndex];
        
        // 如果gt=0的话代表该位点没有变异，则gt改为0/0
        if (gt == "0") {
            outInformations += "\t0/0";
        } else {  // gt不为0的时候
            // 如果位点覆盖度中有0,字段，代表该位点为纯和的突变位点
            if (siteCoverage.find("0,") < siteCoverage.length()) {
                outInformations += "\t1/1";
            } else {  // 如果位点覆盖度没有0,字段，代表该位点为杂合的突变位点
                outInformations += "\t0/1";
            }
        }

        // 将informations最后的字段加上
        for (int i = 1; i < gtVector.size(); i++) {
            outInformations += ":" + gtVector[i];
        }

        outInformations += "\n";

        // 每10Mb写入一次文件
        if (outInformations.size() > bufferSize_) {
            SaveClass.save(outInformations);

            // 清空字符串
            outInformations.clear();
        }
    }

    // 最后写入一次文件
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

    // check file是否排序
    check_vcf_sort(evaluateFileName_);
    
    // 保存分型正确的结果
    string truetxt;
    SAVE trueFile("genotype.true.vcf.gz");

    // 保存分型错误的结果
    string genotypeMisTxt;
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // 保存miscall的结果
    string misCallTxt;
    SAVE misCallFile("miscall.vcf.gz");

    // 没有找到的变异
    string failCallTxt;
    SAVE failCallFile("failcall.vcf.gz");

    // 定义vector
    vector<uint32_t> genotypeLenVec;
    vector<uint32_t> misgenotypeLenVec;
    vector<uint32_t> miscallLenVec;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // 跳过注释行
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }
        
        // 根据FILTER字段进行过滤
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // 如果不是(PASS)则直接跳过
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // 获取基因型信息
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // 根据基因型进行过滤
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
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

        // 获取ref和单倍型的长度
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // 获取变异的长度
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // 添加长度到rocAllMap中
        rocCallMap_[rocNum].push_back(svLength);

        // 进行基因分型评估
        // 检查染色体是否存在，不存在就是miscall
        auto iter1 = trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM);
        if (iter1 == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                    << "Warning: " << INFOSTRUCTTMP.CHROM << " not in true set.\n";
            miscallLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }

        // 记录genotype的结果
        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype
        uint32_t trueRefLen;
        string trueVcfInfo;
        int32_t trueSvLen;

        // 在对应染色体的map中查找真值的refLen qryLenVec
        auto iter2 = iter1->second.find(INFOSTRUCTTMP.POS);
        if (iter2 != iter1->second.end()) {  // 软件找到了，接下来判断分型是否正确
            // true变异信息
            trueRefLen = get<0>(iter2->second);
            const vector<uint32_t>& trueQryLenVec = get<1>(iter2->second);
            trueVcfInfo = get<2>(iter2->second);
            trueSvLen = get<3>(iter2->second);
            vector<int> trueGtVec = get<4>(iter2->second);

            // 如果长度为0，则报错，并退出代码。
            if (trueQryLenVec.size() == 0) {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: trueQryLenVec.size() == 0 -> " 
                    << INFOSTRUCTTMP.CHROM << " " 
                    << INFOSTRUCTTMP.POS << endl;
                exit(1);
            }

            // 遍历等位基因
            int trueQryLenVecIdx = -1;  // 记录真值中用过的单倍型，防止重复评估
            for (size_t i = 0; i < qryLenVec.size(); i++) {
                uint32_t qryLen = qryLenVec[i];
                int callGt = gtVec[i];

                for (size_t j = 0; j < trueQryLenVec.size(); j++) {
                    // 用过的单倍型跳过
                    if (j == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[j];
                    int trueGt = trueGtVec[j];

                    // callGt和trueGt中有一个是0，但另一个不是，则下一个循环，防止SNP判断时候出错
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // 缺失
                    if (refLen >= 50 && qryLen < 50) {
                        if ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // 插入
                        if ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25)
                        {
                            genotypeTrueNum++;
                            trueQryLenVecIdx = j;
                            
                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // 替换
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
                stop:; // 如果找到了，则退出嵌套循环。结束该单倍型的循环，继续下一个单倍型
            }
        }

        // 判断分型的结果
        if (genotypeTrueNum == 0) {  // 找到的变异不在真集中 
            miscallLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        } else {
            // 分型正确
            if (genotypeTrueNum >= gtVec.size()) {
                genotypeLenVec.push_back(trueSvLen);
                truetxt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // 添加长度到rocTrueMap中
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // 分型错误
                misgenotypeLenVec.push_back(trueSvLen);
                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                // 该判断是用recall来计算roc，如果不开启是只用genotype计算roc
                if (rocType_ == "recall") {
                    // 添加长度到rocTrueMap中
                    rocRecallMap_[rocNum].push_back(svLength);
                }
                
            }
            // 删除已经被用过的变异
            trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
        }

        // save result
        if (truetxt.size() > bufferSize_ || 
            genotypeMisTxt.size() > bufferSize_ || 
            misCallTxt.size() > bufferSize_) // 每10Mb写入一次
        {
            trueFile.save(truetxt);  // 分型正确
            genotypeMisFile.save(genotypeMisTxt);  // 分型错误
            misCallFile.save(misCallTxt);  // 不在真集中的变异

            // 清空字符串
            truetxt.clear();
            genotypeMisTxt.clear();
            misCallTxt.clear();
        }
    }

    // save result
    trueFile.save(truetxt);  // 分型正确
    genotypeMisFile.save(genotypeMisTxt);  // 分型错误
    misCallFile.save(misCallTxt);  // 不在真集中的变异

    // 未找到的真实变异
    for (auto it1 : trueVCFStrcuture_.chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (auto it2 : it1.second) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            failCallTxt += get<2>(it2.second) + "\n";

            if (failCallTxt.size() > bufferSize_) {  // 每10Mb写入一次
                failCallFile.save(failCallTxt);  // 未找到的真实变异

                // 清空字符串
                failCallTxt.clear();
            }
        }
    }
    failCallFile.save(failCallTxt);  // 未找到的真实变异

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

    // 长度结果输出
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

    // 释放内存
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

    // 检查vcf文件是否排序
    check_vcf_sort(evaluateFileName_);

    // 保存正确的结果
    string truetxt;
    SAVE trueFile("recall.true.vcf.gz");

    // 保存分型正确的结果
    string true_Gt_txt;
    SAVE trueGtFile("genotype.true.vcf.gz");

    // 保存分型错误的结果
    string genotypeMisTxt;
    SAVE genotypeMisFile("genotype.err.vcf.gz");

    // 保存miscall的结果
    string misCallTxt;
    SAVE misCallFile("miscall.vcf.gz");

    // 定义vector
    vector<uint32_t> genotypeLenVec;
    vector<uint32_t> callLenVec;
    vector<uint32_t> recallLenVec;

    string preChr; // 用于判断是不是新的染色体，是的话再重新构建vector。

    // 记录已经被用过的变异，防止重复评估
    map<string,vector<uint32_t> > selectChrStartMap;

    // 二分查找法左右索引
    int leftIdxTmp = 0;
    int rightIdxTmp = 0;

    // input file stream
    VCFINFOSTRUCT INFOSTRUCTTMP;
    VCFOPEN VCFOPENCLASS(evaluateFileName_);

    while(VCFOPENCLASS.read(INFOSTRUCTTMP)) {
        // 跳过注释行
        if (INFOSTRUCTTMP.line.find("#") != string::npos) {
            continue;
        }

        // 根据FILTER字段进行过滤
        if (INFOSTRUCTTMP.FILTER != "PASS") {  // 如果不是(PASS)则直接跳过
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.FILTER << "' != 'PASS', skip this site -> "
                << INFOSTRUCTTMP.CHROM << "\t" 
                << INFOSTRUCTTMP.POS << endl;
            continue;
        }

        // 获取基因型信息
        vector<int> gtVec = get_gt(
            INFOSTRUCTTMP.lineVec
        );
        string gt = join(gtVec, "/");
        // 根据基因型进行过滤
        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) {  //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
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


        // 获取ref和单倍型的长度
        uint32_t refLen;
        vector<uint32_t> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            INFOSTRUCTTMP.ID, 
            INFOSTRUCTTMP.REF, 
            INFOSTRUCTTMP.ALT, 
            gtVec, 
            "hap"
        );

        // 获取变异的长度
        int32_t svLength = sv_length_select(
            refLen, 
            qryLenVec, 
            gtVec
        );

        // 添加长度到rocAllMap中
        rocCallMap_[rocNum].push_back(svLength);


        // recall
        // 检查染色体是否存在，不存在就是miscall
        if (trueVCFStrcuture_.chrStartLenInfoGtTupMap.find(INFOSTRUCTTMP.CHROM) == trueVCFStrcuture_.chrStartLenInfoGtTupMap.end()) {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Warning: '" << INFOSTRUCTTMP.CHROM << "' is not in the true set.\n";
            callLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
            continue;
        }
    

        if (preChr != INFOSTRUCTTMP.CHROM) {  // 二分查找的左右索引归零，用于加快查询速度
            // 左右索引归零
            leftIdxTmp = 0;
            rightIdxTmp = 0;
        }

        // 记录二分查找法的索引，多个等位基因用同一个refStart
        int indexLeft = -1;
        int indexRight = -1;

        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype

        // 记录recall的  trueRefStart, trueRefLen, trueSvLen、truevcfInfo
        uint32_t trueRefStart;
        uint32_t trueRefLen;
        int32_t trueSvLen;
        string truevcfInfo;

        // 记录真集中被用过的单倍型
        int trueQryLenVecIdx = -1;

        // 遍历多等位基因
        for (size_t i = 0; i < qryLenVec.size(); i++) {
            uint32_t qryLen = qryLenVec[i];
            int callGt = gtVec[i];

            // 二分查找法找200/10bp内的变异索引。
            if (indexLeft == -1 && indexRight == -1) {  // 新的等位基因再查找
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

            // 遍历二分查找法的索引
            for (int j = indexLeft; j <= indexRight; j++) {
                // true变异的信息
                trueRefStart = trueVCFStrcuture_.refStartVecMap[INFOSTRUCTTMP.CHROM][j];

                // 先检查变异有没有被用过，如果被删除了就下一个坐标
                if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(trueRefStart) == trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                    continue;
                }
                
                trueRefLen = get<0>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                const vector<uint32_t>& trueQryLenVec = get<1>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                truevcfInfo = get<2>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                trueSvLen = get<3>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);
                vector<int> trueGtVec = get<4>(trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM][trueRefStart]);

                for (size_t k = 0; k < trueQryLenVec.size(); k++) {  // 一个位点有多个等位基因的时候，trueSeqLenVec存储了各个等位基因的长度，因此遍历它看该位点的变异是否和merge里边的一致，包含0
                    // 用过的单倍型跳过
                    if (k == trueQryLenVecIdx) {
                        continue;
                    }

                    uint32_t trueQryLen = trueQryLenVec[k];
                    int trueGt = trueGtVec[k];

                    // callGt和trueGt中有一个是0，但另一个不是，则下一个循环，防止SNP判断时候出错
                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0)) {
                        continue;;
                    }

                    // 缺失
                    if (refLen >= 50 && qryLen < 50) {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen < 50 && qryLen >= 50) {  // 插入
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 50 && qryLen >= 50) {  // 替换
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=200) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=200) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen == 1 && qryLen == 1) {  // snp
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else if (refLen >= 3 && qryLen <= 2) {  // indel-del
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }  else if (refLen <= 2 && qryLen >= 3) {  // indel-ins
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=10) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=10) &&
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    } else {
                        if ((abs(static_cast<int32_t>(INFOSTRUCTTMP.POS) - static_cast<int32_t>(trueRefStart))<=1) &&
                            (abs((static_cast<int32_t>(INFOSTRUCTTMP.POS) + static_cast<int32_t>(refLen))-(static_cast<int32_t>(trueRefStart) + static_cast<int32_t>(trueRefLen)))<=1) &&
                            ((std::abs(static_cast<int32_t>(refLen) - static_cast<int32_t>(trueRefLen)) / static_cast<float>(trueRefLen)) <= 0.25) && 
                            ((std::abs(static_cast<int32_t>(qryLen) - static_cast<int32_t>(trueQryLen)) / static_cast<float>(trueQryLen)) <= 0.25))
                        {
                            // 记录单倍型recall and genotype 正确
                            genotypeTrueNum++;

                            // 多个等位基因用同一个起始位置
                            indexLeft = j;
                            indexRight = j;

                            // 记录用过单倍型的索引
                            trueQryLenVecIdx = k;

                            goto stop;
                        }
                    }
                }
            }
            stop:; // 如果找到了，则退出嵌套循环。结束该单倍型的循环，继续下一个单倍型
        }

        // 判断寻找的结果并添加 
        if (genotypeTrueNum > 0) {  // 大于0代表找到了
            // 软件call的vcf-vector
            callLenVec.push_back(trueSvLen);

            truetxt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
            recallLenVec.push_back(trueSvLen);

            // 分型正确
            if (genotypeTrueNum >= gtVec.size()) {  // 大于等于gtVec.size()代表分型正确
                genotypeLenVec.push_back(trueSvLen);
                true_Gt_txt += "recall_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                // 添加长度到rocTrueMap中
                rocRecallMap_[rocNum].push_back(svLength);
            } else {  // 分型错误
                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n" + 
                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                // 该判断是用recall来计算roc，如果不开启是只用genotype计算roc
                if (rocType_ == "recall")
                {
                    // 添加长度到rocTrueMap中
                    rocRecallMap_[rocNum].push_back(svLength);
                }
            }

            // 删除已经被用过的变异
            if (trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].find(INFOSTRUCTTMP.POS) != trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].end()) {
                trueVCFStrcuture_.chrStartLenInfoGtTupMap[INFOSTRUCTTMP.CHROM].erase(INFOSTRUCTTMP.POS);
            }
        } else {  // 如果上边循环没找到真集里的变异，则用软件自己找到的长度来添加
            callLenVec.push_back(svLength);
            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + INFOSTRUCTTMP.line + "\n";
        }

        // save result
        if (truetxt.size() > bufferSize_ || 
            true_Gt_txt.size() > bufferSize_ || 
            genotypeMisTxt.size() > bufferSize_ || 
            misCallTxt.size() > bufferSize_) // 每10Mb写入一次
        {
            trueFile.save(truetxt);
            trueGtFile.save(true_Gt_txt);
            genotypeMisFile.save(genotypeMisTxt);  // 分型错误
            misCallFile.save(misCallTxt);  // 不在真集中的变异

            // 清空字符串
            truetxt.clear();
            true_Gt_txt.clear();
            genotypeMisTxt.clear();
            misCallTxt.clear();
        }
    }

    // save result
    trueFile.save(truetxt);
    trueGtFile.save(true_Gt_txt);
    genotypeMisFile.save(genotypeMisTxt);  // 分型错误
    misCallFile.save(misCallTxt);  // 不在真集中的变异

    // 保存没有找到变异，假阴性
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

    // 长度结果输出
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

    // 释放内存
    outFile.close();

    // 计算ROC
    roc_calculate(trueVCFStrcuture_.allLengthList);
}


/**
 * 获取位点基因型列表.
 *
 * @param lineVec          lineVec
 * @param sampleIdx        sample基因型的索引,默认值0代表最后一列
 * 
 * 
 * @return gtVec           vector <int>
**/
vector<int> VCFRecall::get_gt(
    const vector<string> & lineVec, 
    int sampleIdx
) {
    vector <int> gtVec;  // 位点分型的vector

    // FORMAT字符拆分
    vector<string> formatVec = split(lineVec[8], ":");

    int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // 获取GT的索引位置

    if (gtIndex == formatVec.size()) {  // 判断index是否存在，不存在的话返回基因型都为0。
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << lineVec[0] << ":" << lineVec[1] << endl;
        gtVec = {0, 0};
    } else {  // 如果存在，则进行保存
        string gt;  // 存储基因型字段

        if (sampleIdx == 0) {  // 没有指定列数就是最后一列
            gt = split(lineVec[lineVec.size()-1], ":")[gtIndex];  // gt字段
        } else {
            gt = split(lineVec[sampleIdx], ":")[gtIndex];  // gt字段
        }
        
        string splitStr;  // gt中的分隔符
        if (gt.find("/") != string::npos) {  // 判断‘/’分隔符
            splitStr = "/";
        } else if (gt.find("|") != string::npos) {  // 判断‘|’为分隔符
            splitStr = "|";
        } else {  // 不知道的时候为返回空值
            gtVec = {0, 0};
            return gtVec;
        }

        // 找到gt后，对其按splitStr拆分并循环
        auto splitResult = split(gt, splitStr);
        for (auto it : splitResult) {
            // 如果为'.'，跳过该位点
            if (it == ".") {
                gtVec = {0, 0};
                return gtVec;
            }
            // Prevent error in stof(it) conversion
            try {
                gtVec.push_back(stoi(it));  // 添加到vector中
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
 * 获取位点基因型列表.
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

    if (rocColNum_ == 8) {  // FORMAT字段
        vector<string> lastColVec = split(INFOSTRUCTTMP.lineVec.back(), ":");

        vector<string> rocInfVec = split(INFOSTRUCTTMP.FORMAT, ":");
        // rocNum的下标
        vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec_[1]);

        if (rocItera == rocInfVec.end()) {  // 检查该字段中有没有对应的rocNum信息
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: " << rocKeyVec_[1] << " not in " << rocKeyVec_[0] << " column -> " 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        return stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
    } else {  // INFO字段
        smatch patternResult; // 正则表达式的结果
        regex pattern(rocKeyVec_[1] + "=(\\d+)");

        if (!regex_search(INFOSTRUCTTMP.INFO, patternResult, pattern)) {  // 检查是否匹配到 rocKeyVec_[1] 的值
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
 * 获取变异的长度信息
 *
 * @param refLen            ref长度
 * @param qryLenVec         qry长度列表
 * @param gtVec             基因型列表
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
        if (gtVec[i] == 0) {  // 如果基因型是0，跳过
            continue;
        } else {
            svLength = qryLenVec[i] - refLen;
        }
    }
    
    return svLength;
}


/**
 * 统计变异长度信息
 *
 * @param lengthVec_          划分的区间
 * @param length_list        长度列表
 * 
 * 
 * @return vector<uint64_t>  每个区间的长度
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
 * 获取单倍型对应的长度信息.
 *
 * @param svType                         变异类型
 * @param refSeq                         ref列信息
 * @param qrySeqs                        qry列信息
 * @param gtVec                          位点的分型信息
 * @param lenType                        只取单倍型对应的长度还是所有的长度(hap/all)
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
    // 检查模式是否正确，不正确退出代码
    if (lenType != "hap" && lenType != "all") {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "Error: lenType -> " 
            << lenType << endl;
        exit(1);
    }
    
    uint32_t refLen = refSeq.size();  // ref序列长度
    vector<string> qrySeqVec = split(qrySeqs, ",");  // qry的序列列表
    
    // 检查索引是否越界，如果只要hap的长度时候再检查
    if (lenType == "hap") {
        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
        if (maxGtNum > qrySeqVec.size()) {  // 先检查是否数组是否越界
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: number of genotyping and ALT sequences do not match -> " 
                << qrySeqs << endl;
            exit(1);
        }
    }

    // 构造ref和qry共同的长度索引
    vector<uint32_t> seqLenVec;
    seqLenVec.push_back(refLen);

    // 遍历等位基因列表
    for (size_t i = 0; i < qrySeqVec.size(); i++) {
        string qrySeq = qrySeqVec[i];

        // 临时存储长度
        uint32_t refLenTmp = refLen;
        uint32_t qryLenTmp = qrySeq.size(); 

        if (qrySeq.find(">") != string::npos) {  // GraphTyper2的结果
            // 使用正则表达式提取长度信息
            std::regex reg("SVSIZE=(\\d+)");
            std::smatch match;
            // 检查qry序列中有没有长度信息，没有的话跳过该位点
            if (std::regex_search(qrySeq, match, reg)) {  // <INS:SVSIZE=97:BREAKPOINT1>
                // <INS:SVSIZE=90:BREAKPOINT1>
                if (qrySeq.find("<INS") != string::npos && qrySeq.find("<DEL") == string::npos) {  // 字符段中包含插入
                    refLenTmp = refLen;
                    qryLenTmp = std::stoul(match[1].str());;
                } else if (qrySeq.find("<INS") == string::npos && qrySeq.find("<DEL") != string::npos) {  // 字符段中包含缺失: <DEL:SVSIZE=1720:BREAKPOINT>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLen;
                } else if (qrySeq.find("<DUP") != string::npos) {  // 字符段中包含重复的字段（GraphTyper2）<DUP:SVSIZE=2806:BREAKPOINT1>
                    refLenTmp = std::stoul(match[1].str());;
                    qryLenTmp = refLenTmp * 2;
                }
            } else {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: No length information in column 4 -> " << qrySeq << endl;
                refLenTmp = refLen;
                qryLenTmp = qrySeq.size();
            }

            seqLenVec[0] = refLenTmp; // 重置ref的长度
        }

        // 判断BayesTyper的结果为duplication，且ref_len为1-2，qry_seq却很长时
        if (svType.find("Duplication") != string::npos && refLenTmp <= 2) {  // BayesTyper会把duplication变成插入，所以重新计算长度
            refLenTmp = qryLenTmp;
            qryLenTmp *= 2;
        }

        // 添加qry的长度
        if (seqLenVec.size() == (i + 1)) {  // 判断单倍型确实在自己的位置上，
            seqLenVec.push_back(qryLenTmp);  // 添加qry的长度
        } else {  // 索引和列表长度不符合时报错
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Warning: Incorrect index for haplotype length -> " << qrySeqs << endl;
            exit(1);
        }
    }

    // 找基因型对应的qry长度
    vector<uint32_t> qryLenVec;
    if (lenType == "hap") {  // 只取单倍型的序列长度
        if (gtVec.size() == 1) {  // gt只有一个，代表为纯合的变异，两个单倍型添加一样的长度
            if (gtVec[0] == 0) {  // 如果基因型为0，跳过该位点
                // 返回 '0/0'
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
    } else {  // 取所有的长度
        for (size_t i = 1; i < seqLenVec.size(); i++) {
            qryLenVec.push_back(seqLenVec[i]);
        }
    }

    return make_tuple(seqLenVec[0], qryLenVec);
}


/**
 * 保存recall中failCall的结果
 *
 * @param chrStartLenInfoGtTupMap     真集中用剩下的vcf信息
 * @param outFileName                 输出文件名
 * 
 * 
 * @return int             0
**/
int VCFRecall::saveFailCall(
    map<string, map<uint32_t, tuple<uint32_t, vector<uint32_t>, string, int32_t, vector<int> > > > & chrStartLenInfoGtTupMap, 
    const string & outFileName
)
{
    // 没有找到的变异
    string failCallTxt;
    // 输出文件流
    gzFile gzfp = gzopen(outFileName.c_str(), "wb");

    // 遍历字典，将没有找到的vcf进行保存
    for (const auto& [_, StartLenInfoGtTupMap] : chrStartLenInfoGtTupMap) {  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        for (const auto& [_, LenInfoGtTup] : StartLenInfoGtTupMap) {  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            failCallTxt += get<2>(LenInfoGtTup) + "\n";

            // save result
            if (failCallTxt.size() > bufferSize_) {  // 每10Mb写入一次
                gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

                // 清空字符串
                failCallTxt.clear();
            }
        }
    }
    // save result
    gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

    // 释放内存，关闭文件
    gzclose(gzfp);

    return 0;
}


/**
 * 保存recall中failCall的结果
 *
 * @param allLengthList     所有变异长度
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
    
    // 计算所有的长度计算roc
    // 从大往小循环
    string outRoc = rocType_ + "." + rocKey_ + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
    uint64_t allTruePositives = allLengthList.size(); // 所有正确的位点数量
    uint64_t truePositives = 0; // 真阳性
    uint64_t falsePositives = 0; // 假阳性
    uint64_t trueNegatives = 0; // 真阴性
    uint64_t falseNegatives = 0; // 假阴性

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
    allTruePositives = 0; // 所有正确的位点数量
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] == 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // 真阳性
    falsePositives = 0; // 假阳性
    trueNegatives = 0; // 真阴性
    falseNegatives = 0; // 假阴性

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint32_t allPositiveTmp = 0;
        // all中的snp数量
        // 遍历该score下的长度vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] == 0) {
                allPositiveTmp++;
            }
        }

        uint32_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // 遍历该score下的长度vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] == 0) {
                    truePositiveTmp++;
                }
            }
        }
        // 如果个数是0，跳过该score
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
    allTruePositives = 0; // 所有正确的位点数量
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] > -50 && allLengthList[i] < 50 && allLengthList[i] != 0) {
            allTruePositives++;
        }
    }
    truePositives = 0; // 真阳性
    falsePositives = 0; // 假阳性
    trueNegatives = 0; // 真阴性
    falseNegatives = 0; // 假阴性

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // all中的snp数量
        // 遍历该score下的长度vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] > -50 && iter1->second[i] < 50 && iter1->second[i] != 0) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // 遍历该score下的长度vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0) {
                    truePositiveTmp++;
                }
            }
        }
        // 如果个数是0，跳过该score
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
    allTruePositives = 0; // 所有正确的位点数量
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] <= -50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // 真阳性
    falsePositives = 0; // 假阳性
    trueNegatives = 0; // 真阴性
    falseNegatives = 0; // 假阴性

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        int allPositiveTmp = 0;
        // all中的snp数量
        // 遍历该score下的长度vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] <= -50) {
                allPositiveTmp++;
            }
        }

        int truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // 遍历该score下的长度vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] <= -50) {
                    truePositiveTmp++;
                }
            }
        }
        // 如果个数是0，跳过该score
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
    allTruePositives = 0; // 所有正确的位点数量
    for (size_t i = 0; i < allLengthList.size(); i++) {
        if (allLengthList[i] >= 50) {
            allTruePositives++;
        }
    }
    truePositives = 0; // 真阳性
    falsePositives = 0; // 假阳性
    trueNegatives = 0; // 真阴性
    falseNegatives = 0; // 假阴性

    for (auto iter1 = rocCallMap_.rbegin(); iter1 != rocCallMap_.rend(); iter1++) {
        uint64_t allPositiveTmp = 0;
        // all中的snp数量
        // 遍历该score下的长度vector
        for (size_t i = 0; i < iter1->second.size(); i++) {
            if (iter1->second[i] >= 50) {
                allPositiveTmp++;
            }
        }

        uint64_t truePositiveTmp = 0;
        auto findIter1 = rocRecallMap_.find(iter1->first);
        if (findIter1 != rocRecallMap_.end()) {
            // 遍历该score下的长度vector
            for (size_t i = 0; i < findIter1->second.size(); i++) {
                if (findIter1->second[i] >= 50) {
                    truePositiveTmp++;
                }
            }
        }
        // 如果个数是0，跳过该score
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