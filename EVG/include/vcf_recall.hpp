#ifndef vcf_recall_hpp
#define vcf_recall_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <cstdlib>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"

using namespace std;

namespace RECALL
{
    // 全局变量
    // rocKey  
    string rocKey;  // (FORMAT.<name>/INFO.<name>) [FORMAT.DP]
    // rocType判断是用recall来计算roc，还是只用genotype计算roc  recall/genotype [genotype]
    string rocType;
    

    struct vcfStructure
    {
        map<string,vector<int> > refStartVecMap;  // map<chr,vector<refStart> >
        
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > chrStartLenInfoGtVecMap;  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >
        
        vector<int> allLengthList;  // 所有变异的长度
    };


    /**
	 * roc参数修改.
	 *
	 * @param rocKeyTmp        main传入的rocKey
     * @param rocTypeTmp       main传入的rocTypeTmp
     * 
     * 
     * @return int             0
	**/
    int change_roc(
        string rocKeyTmp, 
        string rocTypeTmp
    )
    {
        rocKey = rocKeyTmp;
        rocType = rocTypeTmp;
        return 0;
    }


    // 对变异长度进行统计的Bool函数
    class f_mod
    {
    private:
        int dv1;
        int dv2;

    public:
        f_mod(int d1 = 1, int d2 = 2) : dv1(d1), dv2(d2) {}
    
        bool operator() (int x) {return  dv1 <= x && x <= dv2;}
    };


    int sv_length_select(
        const int & refLen, 
        const vector<int> & qryLenVec, 
        const vector<int> & gtVec
    );
    vcfStructure build_index(
        string genotype_filename
    );
    vector<int> count_num(
        vector<string> sv_length, 
        vector<int> length_list
    );
    string gramtools_convert(
        string evaluateFilename
    );
    int evulate_gt(
        string evaluateFilename, 
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        vector<int> & all_length_list
    );
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx = 0
    );
    tuple<int, vector<int> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );
    int saveFailCall(
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        const string & outFileName
    );
    void roc_calculate(
        const map<float, vector<int> > & rocCallMap, 
        const map<float, vector<int> > & rocRecallMap, 
        const vector<int> & all_length_list
    );


    

    /**
	 * 对vcf文件进行排序并构建索引
	 *
	 * @param genotype_filename    输入文件
     * 
     * 
     * @return outStructure        vcfStructure
	**/
    vcfStructure build_index(
        string genotype_filename
    )
    {
        vcfStructure outStructure;

        // 检查vcf文件是否排序
        check_vcf_sort(genotype_filename);

        // 输入文件流
        gzFile gzfp = gzopen(genotype_filename.c_str(), "rb");

        // 循环将genotype添加到genotype vector中
        vector<int> allLengthList;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Build index.\n";

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "'" << genotype_filename << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // 一次只读1024字节的数据

            while(gzgets(gzfp, line, 1024))
            {
                informations += line;

                if (informations.find("\n") != string::npos) // 一行结束
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // 去掉换行符

                    // 跳过注释行
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo拆分
                        string chromosome = informationsVec[0];  // 染色体信息
                        int refStart = stoi(informationsVec[1]);  // 变异的起始位置
                        string svType = informationsVec[2];  // 判断BayesTyper的结果为duplication
                        string refSeq = informationsVec[3];  // 变异的ref序列
                        string qrySeqs = informationsVec[4];  // 变异的qry序列

                        vector<int> gtVec = get_gt(
                            informationsVec
                        );

                        string gt = join(gtVec, "/");

                        // 根据基因型进行过滤
                        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
                        {
                            // cerr << "[" << __func__ << "::" << getTime() << "] " 
                            //      << "Warning: skip -> " << informations << endl;
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // 获取ref和单倍型的长度
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // 获取变异的长度
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // 保存所有变异的长度
                        outStructure.allLengthList.push_back(svLength);

                        // 变异信息
                        outStructure.refStartVecMap[chromosome].push_back(refStart);

                        // 先初始化哈希表
                        if (outStructure.chrStartLenInfoGtVecMap.find(chromosome) == outStructure.chrStartLenInfoGtVecMap.end())
                        {
                            outStructure.chrStartLenInfoGtVecMap[chromosome];
                        }
                        outStructure.chrStartLenInfoGtVecMap[chromosome][refStart] = make_tuple(refLen, qryLenVec, informations, svLength, gtVec);
                    }

                    // 清空字符串
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // 释放内存，关闭文件
        gzclose(gzfp);

        return outStructure;
    }


    /**
	 * 对gramtools的结果进行转换
	 *
	 * @param evaluateFilename    输入文件
     * 
     * 
     * @return outFileName        输出文件
	**/
    string gramtools_convert(
        string evaluateFilename
    )
    {
        string informations;
        string outInformationsTitle;
        string outInformations;

        // 检查文件是否排序
        check_vcf_sort(
            evaluateFilename
        );

        // 保存结果
        // 输出文件流
        vector<string> prefixVec = split(evaluateFilename, "/");
        string outFileName = "convert." + prefixVec[prefixVec.size()-1];
        if (outFileName.find(".gz") == string::npos && outFileName.find(".GZ") == string::npos)
        {
            outFileName += ".gz";
        }
        
        gzFile gzfp1 = gzopen(outFileName.c_str(), "wb");

        // 打开文件
        if(!gzfp1)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << outFileName 
                << "': No such file or directory." 
                << endl;
            exit(1);
        }


        // 输入文件流
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // 一次只读1024字节的数据

            while(gzgets(gzfp, line, 1024))
            {
                informations += line;

                if (informations.find("\n") != string::npos) // 一行结束
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // 去掉换行符

                    // 跳过注释行
                    if (informations.find("#") == string::npos)
                    {
                        // 保存注释行
                        if (outInformationsTitle.size() > 0)
                        {
                            gzwrite(gzfp1, outInformationsTitle.c_str(), outInformationsTitle.length());

                            // 清空注释行字符串
                            outInformationsTitle.clear();
                            string().swap(outInformationsTitle);
                        }
                        // 每10Mb写入一次文件
                        if (outInformations.size() > 10000000)
                        {
                            gzwrite(gzfp1, outInformations.c_str(), outInformations.length());

                            // 清空字符串
                            outInformations.clear();
                            string().swap(outInformations);
                        }

                        // 对字符串进行分割并保存到vector中
                        vector<string> informationsVec = split(informations, "\t");

                        // 提取informationsVector中基因型和位点覆盖度信息
                        vector<string> gtVector = split(strip(informationsVec[informationsVec.size()-1], '\n'), ":");
                        vector<string> formatVector = split(strip(informationsVec[8], '\n'), ":");
                        string gt = gtVector[0];

                        // 找FORMAT字段中COV的位置
                        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
                        int covIndex = 0;
                        if (covItera != formatVector.end()) // FORMAT中有COV
                        {
                            covIndex = distance(formatVector.begin(), covItera);
                        }
                        else // 没有的话直接转换，继续下一个循环
                        {
                            outInformations += informations + "\n";
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // FILTER字段前边的内容不变，直接加到outInformations上
                        for (int i = 0; i < informationsVec.size()-1; i++)
                        {
                            if (i == 6) // FILTER列改为PASS，源文件有的为.
                            {
                                outInformations += "PASS\t";
                            }
                            else
                            {
                                outInformations += informationsVec[i] + "\t";
                            }
                        }

                        string siteCoverage = gtVector[covIndex];
                        

                        // 如果gt=0的话代表该位点没有变异，则gt改为0/0
                        if (gt == "0")
                        {
                            outInformations += "\t0/0";
                        }
                        // gt不为0的时候
                        else
                        {
                            // 如果位点覆盖度中有0,字段，代表该位点为纯和的突变位点
                            if (siteCoverage.find("0,") < siteCoverage.length())
                            {
                                outInformations += "\t1/1";
                            }
                            // 如果位点覆盖度没有0,字段，代表该位点为杂合的突变位点
                            else
                            {
                                outInformations += "\t0/1";
                            }
                        }

                        // 将informations最后的字段加上
                        for (int i = 1; i < gtVector.size(); i++)
                        {
                            outInformations += ":" + gtVector[i];
                        }

                        outInformations += "\n";
                    }
                    // 保存注释行信息
                    else
                    {
                        outInformationsTitle += informations + "\n";
                    }
                    
                    // 清空字符串
                    informations.clear();
                    string().swap(informations);
                }
            }
        }

        // 释放内存，关闭文件
        gzclose(gzfp);

        // 最后写入一次文件
        gzwrite(gzfp1, outInformations.c_str(), outInformations.length());
        // 释放内存，关闭文件
        gzclose(gzfp1);
        
        return outFileName;
    }


    /**
	 * genotype评估函数
	 *
	 * @param evaluateFilename            输入文件（待评估）
     * @param chrStartLenInfoGtVecMap     真集所有的vcf信息
     * @param all_length_list             真集所有vcf长度
     * 
     * 
     * @return int             0
	**/
    int evulate_gt(
        string evaluateFilename, 
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        vector<int> & all_length_list
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Genotype." << endl;
        
        // 变异长度范围
        vector<string> sv_length = {
            "-999999/-10000", 
            "-9999/-5000", 
            "-4999/-2500", 
            "-2499/-1000",
            "-999/-500", 
            "-499/-400", 
            "-399/-300", 
            "-299/-200", 
            "-199/-100", 
            "-99/-50", 
            "-49/-1", 
            "0/0", 
            "1/49", 
            "50/99", 
            "100/199", 
            "200/299", 
            "300/399", 
            "400/499", 
            "500/999", 
            "1000/2499",  
            "2500/4999", 
            "5000/9999",
            "10000/999999"
        };

        // ROC -> save the number of DP or GQ 
        map<float, vector<int> > rocCallMap; // map<DP/GQ, vector<length>> 软件找到所有的
        map<float, vector<int> > rocRecallMap; // map<DP/GQ, vector<length>> 软件找到正确的
        vector<string> rocKeyVec = split(rocKey, "."); // vector<colname, key>
        int rocColNum = 0;
        if (rocKeyVec[0] == "INFO")
        {
            rocColNum = 7;
        }
        else if (rocKeyVec[0] == "FORMAT")
        {
            rocColNum = 8;
        }
        else
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: please enter the correct column name to calculate the roc value (INFO/FORMAT): -> " <<  rocKeyVec[0]
                << endl;
            exit(1);
        }

        // 检查文件是否排序
        check_vcf_sort(evaluateFilename);
        
        // 保存分型正确的结果
        string true_txt;
        gzFile trueFile = gzopen("genotype.true.vcf.gz", "wb");
        // 打开文件
        if(!trueFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 保存分型错误的结果
        string genotypeMisTxt;
        gzFile genotypeMisFile = gzopen("genotype.err.vcf.gz", "wb");
        // 打开文件
        if(!genotypeMisFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.err.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 保存miscall的结果
        string misCallTxt;
        gzFile misCallFile = gzopen("miscall.vcf.gz", "wb");
        // 打开文件
        if(!misCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'miscall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 没有找到的变异
        string failCallTxt;
        gzFile failCallFile = gzopen("failcall.vcf.gz", "wb");
        // 打开文件
        if(!failCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'failcall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 定义vector
        vector<int> genotype_length_list;
        vector<int> misgenotype_length_list;
        vector<int> miscall_length_list;

        // 输入文件流
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // 一次只读1024字节的数据

            while(gzgets(gzfp, line, 1024))
            {
                informations += line;

                if (informations.find("\n") != string::npos) // 一行结束
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // 去掉换行符

                    // 跳过注释行
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo拆分
                        string chromosome = informationsVec[0];
                        int refStart = stoi(informationsVec[1]);  // 变异的起始位置
                        string svType = informationsVec[2];  // 判断BayesTyper的结果为duplication
                        string refSeq = informationsVec[3];  // 变异的ref序列
                        string qrySeqs = informationsVec[4];  // 变异的qry序列
                        string filter = informationsVec[6];

                        // 根据FILTER字段进行过滤
                        if (filter != "PASS") // 如果不是(PASS)则直接跳过
                        {
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // 获取基因型信息
                        vector<int> gtVec = get_gt(
                            informationsVec
                        );
                        string gt = join(gtVec, "/");
                        // 根据基因型进行过滤
                        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
                        {
                            // cerr << "[" << __func__ << "::" << getTime() << "] " 
                            //      << "Warning: skip -> " << informations << endl;
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        /* --------------------------------- roc ------------------------------------ */
                        vector<string> lastColVec = split(informationsVec[informationsVec.size()-1], ":");
                        float rocNum;
                        if (rocColNum == 8) // FORMAT字段
                        {
                            vector<string> rocInfVec = split(strip(informationsVec[rocColNum], '\n'), ":");
                            // rocNum的下标
                            vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec[1]);
                            if (rocItera == rocInfVec.end()) // 检查该字段中有没有对应的rocNum信息
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
                        }
                        else // INFO字段
                        {
                            smatch patternResult; // 正则表达式的结果
                            string rocInfString = strip(informationsVec[rocColNum], '\n');
                            regex pattern(rocKeyVec[1]  + "=(\\d+)");
                            string::const_iterator iterStart = rocInfString.begin();
                            string::const_iterator iterEnd = rocInfString.end();
                            regex_search(iterStart, iterEnd, patternResult, pattern);
                            string rocNumString = patternResult[1];

                            // 检查有没有结果，没有的话报错
                            if (rocNumString.empty())
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(rocNumString);
                        }
                        /* --------------------------------- roc ------------------------------------ */

                        // 获取ref和单倍型的长度
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // 获取变异的长度
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // 添加长度到rocAllMap中
                        rocCallMap[rocNum].push_back(svLength);

                        // 进行基因分型评估
                        // 检查染色体是否存在，不存在就是miscall
                        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > >::iterator iter1 = chrStartLenInfoGtVecMap.find(chromosome);
                        if (iter1 == chrStartLenInfoGtVecMap.end())
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Warning: " << chromosome << " not in chrStartLenInfoGtVecMap.\n";
                            miscall_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // 记录genotype的结果
                        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype
                        int trueRefLen;
                        string trueVcfInfo;
                        int trueSvLen;

                        // 在对应染色体的map中查找真值的refLen qryLenVec
                        map<int, tuple<int, vector<int>, string, int, vector<int> > >::iterator iter2 = iter1->second.find(refStart);
                        if (iter2 != iter1->second.end())  // 软件找到了，接下来判断分型是否正确
                        {
                            // true变异信息
                            trueRefLen = get<0>(iter2->second);
                            vector<int> trueQryLenVec = get<1>(iter2->second);
                            trueVcfInfo = get<2>(iter2->second);
                            trueSvLen = get<3>(iter2->second);
                            vector<int> trueGtVec = get<4>(iter2->second);

                            // 如果长度为0，则报错，并退出代码。
                            if (trueQryLenVec.size() == 0)
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " 
                                    << "Error: trueQryLenVec.size() == 0 -> chromosome:" << chromosome 
                                    << " refStart:" << refStart << endl;
                                exit(1);
                            }

                            // 遍历等位基因
                            int trueQryLenVecIdx = -1;  // 记录真值中用过的单倍型，防止重复评估
                            for (size_t i = 0; i < qryLenVec.size(); i++)
                            {
                                int qryLen = qryLenVec[i];
                                int callGt = gtVec[i];

                                for (size_t j = 0; j < trueQryLenVec.size(); j++)
                                {
                                    // 用过的单倍型跳过
                                    if (j == trueQryLenVecIdx)
                                    {
                                        continue;
                                    }

                                    int trueQryLen = trueQryLenVec[j];
                                    int trueGt = trueGtVec[j];

                                    // callGt和trueGt中有一个是0，但另一个不是，则下一个循环，防止SNP判断时候出错
                                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0))
                                    {
                                        continue;;
                                    }

                                    // 缺失
                                    if (refLen >= 50 && qryLen < 50)
                                    {
                                        if ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;

                                            goto stop;
                                        }
                                    }
                                    // 插入
                                    else if (refLen < 50 && qryLen >= 50)
                                    {
                                        if ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    // 替换
                                    else if (refLen >= 50 && qryLen >= 50)
                                    {
                                        if (((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    // snp
                                    else if (refLen == 1 && qryLen == 1)
                                    {
                                        if (((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    // del
                                    else if (refLen >= 3 && qryLen <= 2)
                                    {
                                        if ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    // ins
                                    else if (refLen <= 2 && qryLen >= 3)
                                    {
                                        if ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    else
                                    {
                                        if ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25 && 
                                        (abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25)
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
                        if (genotypeTrueNum == 0)  // 找到的变异不在真集中
                        {
                            miscall_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                        else
                        {
                            // 分型正确
                            if (genotypeTrueNum >= gtVec.size())
                            {
                                genotype_length_list.push_back(trueSvLen);
                                true_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                                // 添加长度到rocTrueMap中
                                rocRecallMap[rocNum].push_back(svLength);
                            }
                            // 分型错误
                            else
                            {
                                misgenotype_length_list.push_back(trueSvLen);
                                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                                // 该判断是用recall来计算roc，如果不开启是只用genotype计算roc
                                if (rocType == "recall")
                                {
                                    // 添加长度到rocTrueMap中
                                    rocRecallMap[rocNum].push_back(svLength);
                                }
                                
                            }
                            // 删除已经被用过的变异
                            chrStartLenInfoGtVecMap[chromosome].erase(refStart);
                        }
                    }
                    // 清空字符串
                    informations.clear();
                    string().swap(informations);

                    // 保存结果
                    if (true_txt.size() > 10000000 || 
                        genotypeMisTxt.size() > 10000000 || 
                        misCallTxt.size() > 10000000) // 每10Mb写入一次
                    {
                        gzwrite(trueFile, true_txt.c_str(), true_txt.length()); // 分型正确
                        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // 分型错误
                        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // 不在真集中的变异

                        // 清空字符串
                        true_txt.clear();
                        genotypeMisTxt.clear();
                        misCallTxt.clear();
                        string().swap(true_txt);
                        string().swap(genotypeMisTxt);
                        string().swap(misCallTxt);
                    }
                }
            }
        }
        // 释放内存，关闭文件
        gzclose(gzfp);

        // 保存结果
        gzwrite(trueFile, true_txt.c_str(), true_txt.length()); // 分型正确

        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // 分型错误

        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // 不在真集中的变异

        // 未找到的真实变异
        for (auto it1 : chrStartLenInfoGtVecMap)  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        {
            for (auto it2 : it1.second)  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            {
                failCallTxt += get<3>(it2.second) + "\n";

                if (failCallTxt.size() > 10000000) // 每10Mb写入一次
                {
                    gzwrite(failCallFile, failCallTxt.c_str(), failCallTxt.length()); // 分型正确

                    // 清空字符串
                    failCallTxt.clear();
                    string().swap(failCallTxt);
                }
            }
        }
        gzwrite(failCallFile, failCallTxt.c_str(), failCallTxt.length());

        // 关闭文件
        gzclose(trueFile);
        gzclose(genotypeMisFile);
        gzclose(misCallFile);
        gzclose(failCallFile);

        int sv_genotype_recall = genotype_length_list.size();
        int sv_misgenotype_recall = misgenotype_length_list.size();
        int sv_mis_call = miscall_length_list.size();
        int sv_all = all_length_list.size();

        // 保存结果
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

        vector<int> all_length_count;
        vector<int> mis_call_length_count;
        vector<int> genotype_length_count;
        vector<int> misgenotype_length_count;

        all_length_count = count_num(sv_length, all_length_list);
        mis_call_length_count = count_num(sv_length, miscall_length_list);
        genotype_length_count = count_num(sv_length, genotype_length_list);
        misgenotype_length_count = count_num(sv_length, misgenotype_length_list);
        
        outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

        // 长度结果输出
        for (int i = 0; i < sv_length.size(); i++)
        {
            outFile << sv_length[i] << ": " 
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

        roc_calculate(rocCallMap, rocRecallMap, all_length_list);

        return 0;
    }


    /**
	 * recall评估函数
	 *
	 * @param evaluateFilename            输入文件（待评估）
     * @param trueVcfStructure            build_index()输出结果
     * 
     * 
     * @return int             0
	**/
    int evulate_recall(
        const string & evaluateFilename, 
        vcfStructure & trueVcfStructure
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Recall." << endl;

        // 变异长度范围
        vector<string> sv_length = {
            "-999999/-10000", 
            "-9999/-5000", 
            "-4999/-2500", 
            "-2499/-1000",
            "-999/-500", 
            "-499/-400", 
            "-399/-300", 
            "-299/-200", 
            "-199/-100", 
            "-99/-50",
             "-49/-1", 
             "0/0", 
             "1/49", 
             "50/99", 
             "100/199", 
             "200/299", 
             "300/399", 
             "400/499", 
             "500/999", 
             "1000/2499", 
             "2500/4999", 
             "5000/9999", 
             "10000/999999"
        };

        // ROC -> save the number of DP or GQ 
        map<float, vector<int>> rocCallMap; // map<DP/GQ, vector<length>> 软件找到所有的
        map<float, vector<int>> rocRecallMap; // map<DP/GQ, vector<length>> 软件找到正确的
        vector<string> rocKeyVec = split(rocKey, "."); // vector<colname, key>
        int rocColNum = 0;
        if (rocKeyVec[0] == "INFO")
        {
            rocColNum = 7;
        }
        else if (rocKeyVec[0] == "FORMAT")
        {
            rocColNum = 8;
        }
        else
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: please enter the correct column name to calculate the roc value (INFO/FORMAT) -> " <<  rocKeyVec[0]
                << endl;
            exit(1);
        }

        // 检查vcf文件是否排序
        check_vcf_sort(evaluateFilename);

        // 保存正确的结果
        string true_txt{};
        gzFile trueFile = gzopen("recall.true.vcf.gz", "wb");
        // 打开文件
        if(!trueFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'recall.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 保存分型正确的结果
        string true_Gt_txt{};
        gzFile trueGtFile = gzopen("genotype.true.vcf.gz", "wb");
        // 打开文件
        if(!trueGtFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 保存分型错误的结果
        string genotypeMisTxt{};
        gzFile genotypeMisFile = gzopen("genotype.err.vcf.gz", "wb");
        // 打开文件
        if(!genotypeMisFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.err.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 保存miscall的结果
        string misCallTxt{};
        gzFile misCallFile = gzopen("miscall.vcf.gz", "wb");
        // 打开文件
        if(!misCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'miscall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // 定义vector
        vector<int> genotype_length_list;
        vector<int> call_length_list;
        vector<int> recall_length_list;

        // 跳过注释行
        string start_chromosome; // 用于判断是不是新的染色体，是的话再重新构建vector。

        // 记录已经被用过的变异，防止重复评估
        map<string,vector<int>> selectChrStartMap;

        vector<int> start_vector;
        vector<int> ref_len_vector;
        vector<int> qry_len_vector;
        vector<string> gt_vector;
        vector<string> vcf_inf_vector;

        // 二分查找法左右索引
        int leftIdxTmp = 0;
        int rightIdxTmp = 0;

        // 输入文件流
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // 一次只读1024字节的数据

            while(gzgets(gzfp, line, 1024))
            {
                informations += line;

                if (informations.find("\n") != string::npos) // 一行结束
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // 去掉换行符

                    // 跳过注释行
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");

                        string chromosome = informationsVec[0];
                        int refStart = stoi(informationsVec[1]);  // 变异的起始位置
                        string svType = informationsVec[2];  // 判断BayesTyper的结果为duplication
                        string refSeq = informationsVec[3];  // 变异的ref序列
                        string qrySeqs = informationsVec[4];  // 变异的qry序列
                        string filter = informationsVec[6];

                        // 根据FILTER字段进行过滤
                        if (filter != "PASS") // 如果不是(PASS)则直接跳过
                        {
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // 获取基因型信息
                        vector<int> gtVec = get_gt(
                            informationsVec
                        );
                        string gt = join(gtVec, "/");
                        // 根据基因型进行过滤
                        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) //如果是(0/0, .)格式的则直接跳过，或者返回了空列表，不构建索引
                        {
                            // cerr << "[" << __func__ << "::" << getTime() << "] " 
                            //      << "Warning: skip -> " << informations << endl;
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }


                        /* --------------------------------- roc ------------------------------------ */
                        vector<string> lastColVec = split(informationsVec[informationsVec.size()-1], ":");
                        float rocNum;
                        if (rocColNum == 8) // FORMAT字段
                        {
                            vector<string> rocInfVec = split(strip(informationsVec[rocColNum], '\n'), ":");
                            // rocNum的下标
                            vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec[1]);
                            if (rocItera == rocInfVec.end()) // 检查该字段中有没有对应的rocNum信息
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
                        }
                        else // INFO字段
                        {
                            smatch patternResult; // 正则表达式的结果
                            string rocInfString = strip(informationsVec[rocColNum], '\n');
                            regex pattern(rocKeyVec[1]  + "=(\\d+)");
                            string::const_iterator iterStart = rocInfString.begin();
                            string::const_iterator iterEnd = rocInfString.end();
                            regex_search(iterStart, iterEnd, patternResult, pattern);
                            string rocNumString = patternResult[1];

                            // 检查有没有结果，没有的话报错
                            if (rocNumString.empty())
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(rocNumString);
                        }
                        /* --------------------------------- roc ------------------------------------ */


                        // 获取ref和单倍型的长度
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // 获取变异的长度
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // 添加长度到rocAllMap中
                        rocCallMap[rocNum].push_back(svLength);


                        // recall
                        // 检查染色体是否存在，不存在就是miscall
                        if (trueVcfStructure.chrStartLenInfoGtVecMap.find(chromosome) == trueVcfStructure.chrStartLenInfoGtVecMap.end())
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Warning: " << chromosome << " not in chrStartLenInfoGtVecMap.\n";
                            call_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                    

                        if (start_chromosome != chromosome)  // 二分查找的左右索引归零，用于加快查询速度
                        {
                            // 左右索引归零
                            leftIdxTmp = 0;
                            rightIdxTmp = 0;
                        }

                        // 记录二分查找法的索引，多个等位基因用同一个refStart
                        int indexLeft = -1;
                        int indexRight = -1;

                        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype

                        // 记录recall的  trueRefStart, trueRefLen, trueSvLen、truevcfInfo
                        int trueRefStart;
                        int trueRefLen;
                        int trueSvLen;
                        string truevcfInfo;

                        // 记录真集中被用过的单倍型
                        int trueQryLenVecIdx = -1;

                        // 遍历多等位基因
                        for (size_t i = 0; i < qryLenVec.size(); i++)
                        {
                            int qryLen = qryLenVec[i];
                            int callGt = gtVec[i];

                            // 二分查找法找200/10bp内的变异索引。
                            if (indexLeft == -1 && indexRight == -1)  // 新的等位基因再查找
                            {
                                if (refLen <= 49 && qryLen <= 49)
                                {
                                    indexLeft = search_Binary_left(trueVcfStructure.refStartVecMap[chromosome], (refStart-10), leftIdxTmp);
                                    leftIdxTmp = indexLeft;
                                    indexRight = search_Binary_right(trueVcfStructure.refStartVecMap[chromosome], (refStart+10), rightIdxTmp);
                                    rightIdxTmp = indexRight;
                                }
                                else
                                {
                                    indexLeft = search_Binary_left(trueVcfStructure.refStartVecMap[chromosome], (refStart-200), leftIdxTmp);
                                    leftIdxTmp = indexLeft;
                                    indexRight = search_Binary_right(trueVcfStructure.refStartVecMap[chromosome], (refStart+200), rightIdxTmp);
                                    rightIdxTmp = indexRight;
                                }
                                if (indexLeft < 0 || indexRight >= trueVcfStructure.refStartVecMap[chromosome].size())
                                {
                                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: out of index, please check the data or code.\n";
                                    exit(1);
                                }
                            }

                            // 遍历二分查找法的索引
                            for (int j = indexLeft; j <= indexRight; j++)
                            {
                                // true变异的信息
                                trueRefStart = trueVcfStructure.refStartVecMap[chromosome][j];

                                // 先检查变异有没有被用过，如果被删除了就下一个坐标
                                if (trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].find(trueRefStart) == trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].end())
                                {
                                    continue;
                                }
                                
                                trueRefLen = get<0>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                vector<int> trueQryLenVec = get<1>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                truevcfInfo = get<2>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                trueSvLen = get<3>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                vector<int> trueGtVec = get<4>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);

                                for (size_t k = 0; k < trueQryLenVec.size(); k++)  // 一个位点有多个等位基因的时候，trueSeqLenVec存储了各个等位基因的长度，因此遍历它看该位点的变异是否和merge里边的一致，包含0
                                {
                                    // 用过的单倍型跳过
                                    if (k == trueQryLenVecIdx)
                                    {
                                        continue;
                                    }

                                    int trueQryLen = trueQryLenVec[k];
                                    int trueGt = trueGtVec[k];

                                    // callGt和trueGt中有一个是0，但另一个不是，则下一个循环，防止SNP判断时候出错
                                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0))
                                    {
                                        continue;;
                                    }

                                    // 缺失
                                    if (refLen >= 50 && qryLen < 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
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
                                    // 插入
                                    else if (refLen < 50 && qryLen >= 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // 替换
                                    else if (refLen >= 50 && qryLen >= 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // snp
                                    else if (refLen == 1 && qryLen == 1)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // indel-del
                                    else if (refLen >= 3 && qryLen <= 2)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
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
                                    // indel-ins
                                    else if (refLen <= 2 && qryLen >= 3)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    else
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                        if (genotypeTrueNum > 0)  // 大于0代表找到了
                        {
                            // 软件call的vcf-vector
                            call_length_list.push_back(trueSvLen);

                            true_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                            recall_length_list.push_back(trueSvLen);

                            // 分型正确
                            if (genotypeTrueNum >= gtVec.size())  // 大于等于gtVec.size()代表分型正确
                            {
                                genotype_length_list.push_back(trueSvLen);
                                true_Gt_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                                // 添加长度到rocTrueMap中
                                rocRecallMap[rocNum].push_back(svLength);
                            }
                            // 分型错误
                            else
                            {
                                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                                // 该判断是用recall来计算roc，如果不开启是只用genotype计算roc
                                if (rocType == "recall")
                                {
                                    // 添加长度到rocTrueMap中
                                    rocRecallMap[rocNum].push_back(svLength);
                                }
                            }

                            // 删除已经被用过的变异
                            if (trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].find(refStart) != trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].end())
                            {
                                trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].erase(refStart);
                            }
                        }
                        else  // 如果上边循环没找到真集里的变异，则用软件自己找到的长度来添加
                        {
                            call_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                        }
                    }
                    // 清空字符串
                    informations.clear();
                    string().swap(informations);

                    // 保存结果
                    if (true_txt.size() > 10000000 || 
                        true_Gt_txt.size() > 10000000 || 
                        genotypeMisTxt.size() > 10000000 || 
                        misCallTxt.size() > 10000000) // 每10Mb写入一次
                    {
                        gzwrite(trueFile, true_txt.c_str(), true_txt.length());
                        gzwrite(trueGtFile, true_Gt_txt.c_str(), true_Gt_txt.length());
                        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // 分型错误
                        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // 不在真集中的变异

                        // 清空字符串
                        true_txt.clear();
                        true_Gt_txt.clear();
                        genotypeMisTxt.clear();
                        misCallTxt.clear();
                        string().swap(true_txt);
                        string().swap(true_Gt_txt);
                        string().swap(genotypeMisTxt);
                        string().swap(misCallTxt);
                    }
            }
            }
        }
        // 释放内存，关闭文件
        gzclose(gzfp);

        // 保存结果
        gzwrite(trueFile, true_txt.c_str(), true_txt.length());

        gzwrite(trueGtFile, true_Gt_txt.c_str(), true_Gt_txt.length());

        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // 分型错误

        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // 不在真集中的变异

        // 保存没有找到变异，假阴性
        saveFailCall(
            trueVcfStructure.chrStartLenInfoGtVecMap, 
            "failcall.vcf.gz"
        );

        // 关闭文件
        gzclose(trueFile);
        gzclose(trueGtFile);
        gzclose(genotypeMisFile);
        gzclose(misCallFile);

        int all_num = trueVcfStructure.allLengthList.size();
        int call_num = call_length_list.size();
        int recall_num = recall_length_list.size();

        int genotype_recall_number = genotype_length_list.size();

        // 保存结果
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

        vector<int> all_num_count = count_num(sv_length, trueVcfStructure.allLengthList);
        vector<int> call_num_count = count_num(sv_length, call_length_list);
        vector<int> recall_num_count = count_num(sv_length, recall_length_list);

        vector<int> genotype_recall_number_count = count_num(sv_length, genotype_length_list);
        
        outFile << "length/length: genotype_recall/misgenotype_call/recall/mis_call/call/fail_call/all\n";

        // 长度结果输出
        for (int i = 0; i < sv_length.size(); i++)
        {
            outFile << sv_length[i] << ": " 
                    << genotype_recall_number_count[i] << "/" 
                    << recall_num_count[i] - genotype_recall_number_count[i] << "/" 
                    << recall_num_count[i] << "/" 
                    << call_num_count[i] - recall_num_count[i] << "/"
                    << call_num_count[i] << "/"
                    << all_num_count[i] - recall_num_count[i] << "/"
                    << all_num_count[i] << endl;
            
            if ( 10 <= i && i <= 12)
            {
                all_num -= all_num_count[i];
                call_num -= call_num_count[i];
                recall_num -= recall_num_count[i];
                genotype_recall_number -= genotype_recall_number_count[i];
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
        roc_calculate(rocCallMap, rocRecallMap, trueVcfStructure.allLengthList);

        return 0;
    }


    /**
	 * 获取位点基因型列表.
	 *
	 * @param informationsVec  vcfInfoList
     * @param sampleIdx        sample基因型的索引,默认值0代表最后一列
     * 
     * 
     * @return gtVec           vector <int>
	**/
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx
    )
    {
        vector <int> gtVec;  // 位点分型的vector

        vector<string> formatVec;  // FORMAT字符拆分
        int formatIndex = 8; // FORMAT所在列
        formatVec = split(informationsVec[formatIndex], ":");

        int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // 获取GT的索引位置

        if (gtIndex == formatVec.size())  // 判断index是否存在，不存在的话返回基因型都为0。
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
            gtVec = {0, 0};
        }
        else  // 如果存在，则进行保存
        {
            string gt;  // 存储基因型字段

            if (sampleIdx == 0) // 没有指定列数就是最后一列
            {
                gt = split(informationsVec[informationsVec.size()-1], ":")[gtIndex];  // gt字段
            }
            else
            {
                gt = split(informationsVec[sampleIdx], ":")[gtIndex];  // gt字段
            }
            
            string splitStr;  // gt中的分隔符
            if (gt.find("/") != string::npos)  // 判断‘/’分隔符
            {
                splitStr = "/";
            }
            else if (gt.find("|") != string::npos)  // 判断‘|’为分隔符
            {
                splitStr = "|";
            }
            else  // 不知道的时候为返回空值
            {
                gtVec = {0, 0};
                return gtVec;
            }
            
            for (auto it : split(gt, splitStr))  // 找到gt后，对其按splitStr拆分并循环
            {
                if (it == ".")  // 如果为'.'，跳过该位点
                {
                    gtVec = {0, 0};
                    return gtVec;
                }
                gtVec.push_back(stoi(it));  // 添加到vector中
            }
        }

        return gtVec;
    }



    /**
	 * 获取变异的长度信息
	 *
     * @param refLen            ref长度
	 * @param qryLenVec         qry长度列表
     * @param gtVec             基因型列表
     * 
     * 
     * @return int              svLength
	**/
    int sv_length_select(
        const int & refLen, 
        const vector<int> & qryLenVec, 
        const vector<int> & gtVec
    )
    {
        int svLength;

        for (size_t i = 0; i < gtVec.size(); i++)
        {
            if (gtVec[i] == 0)  // 如果基因型是0，跳过
            {
                continue;
            }
            else
            {
                svLength = qryLenVec[i] - refLen;
            }
        }
        
        return svLength;
    }


    /**
	 * 统计变异长度信息
	 *
     * @param sv_length         划分的区间
	 * @param length_list       长度列表
     * 
     * 
     * @return vector<int>      每个区间的长度
	**/
    vector<int> count_num(
        vector<string> sv_length, 
        vector<int> length_list
    )
    {
        vector<int> out_length;
        for (int i = 0; i < sv_length.size(); i++)
        {
            int x1 = std::stoi(split(sv_length[i], "/")[0]);
            int x2 = std::stoi(split(sv_length[i], "/")[1]);

            vector<int>::size_type result = count_if(length_list.begin(), length_list.end(), f_mod(x1, x2));

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
     * @return tuple<int, vector<int> >      tuple<refLen, vector<qryLen> >
	**/
    tuple<int, vector<int> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    )
    {
        // 检查模式是否正确，不正确退出代码
        if (lenType != "hap" && lenType != "all")
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: lenType -> " 
                << lenType << endl;
            exit(1);
        }
        
        int refLen = refSeq.size();  // ref序列长度
        vector<string> qrySeqVec = split(qrySeqs, ",");  // qry的序列列表
        
        // 检查索引是否越界，如果只要hap的长度时候再检查
        if (lenType == "hap")
        {
            int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
            if (maxGtNum > qrySeqVec.size())  // 先检查是否数组是否越界
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: number of genotyping and query sequences do not match -> " 
                    << qrySeqs << endl;
                exit(1);
            }
        }

        // 构造ref和qry共同的长度索引
        vector<int> seqLenVec;
        seqLenVec.push_back(refLen);

        // 遍历等位基因列表
        for (size_t i = 0; i < qrySeqVec.size(); i++)
        {
            string qrySeq = qrySeqVec[i];

            // 临时存储长度
            int refLenTmp = refLen;
            int qryLenTmp = qrySeq.size();

            string::size_type idxSvsize = qrySeq.find(">");
            string::size_type idxIns = qrySeq.find("<INS");  // <INS:SVSIZE=90:BREAKPOINT1>
            string::size_type idxDel = qrySeq.find("<DEL");  // <DEL:SVSIZE=1720:BREAKPOINT>
            string::size_type idxDup = qrySeq.find("<DUP");  // <DUP:SVSIZE=10001:COVERAGE>

            // 判断BayesTyper的结果为duplication
            string::size_type idxSvTypeDup = svType.find("Duplication");

            if (idxSvsize != string::npos)  // GraphTyper2的结果
            {
                if (idxIns != string::npos && idxDel == string::npos) // 字符段中包含插入
                {
                    refLenTmp = refLen;
                    qryLenTmp = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                }
                else if (idxIns == string::npos && idxDel != string::npos) // 字符段中包含缺失
                {
                    refLenTmp = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    qryLenTmp = refLen;
                }
                else if (idxDup != string::npos) // 字符段中包含重复的字段（GraphTyper2）<DUP:SVSIZE=2806:BREAKPOINT1>
                {
                    refLenTmp = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    qryLenTmp = refLenTmp * 2;
                }
                else
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: unknown string -> ref_seq:" << refSeq << " qey_seq:" << qrySeq << endl;
                    refLenTmp = refLen;
                    qryLenTmp = stoi(split(split(qrySeqs, "=")[1], ":")[0]);
                }

                seqLenVec[0] = refLenTmp; // 重置ref的长度
            }

            // 判断BayesTyper的结果为duplication，且ref_len为1-2，qry_seq却很长时
            if (idxSvTypeDup != string::npos && refLenTmp <= 2) // BayesTyper会把duplication变成插入，所以重新计算长度
            {
                refLenTmp = qryLenTmp;
                qryLenTmp *= 2;
            }

            // 添加qry的长度
            if (seqLenVec.size() == (i + 1))  // 判断单倍型确实在自己的位置上，
            {
                seqLenVec.push_back(qryLenTmp);  // 添加qry的长度
            }
            else  // 索引和列表长度不符合时报错
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: wrong index for haplotype length -> " 
                    << qrySeqs 
                    << endl;
                exit(1);
            }
        }

        // 找基因型对应的qry长度
        vector<int> qryLenVec;
        if (lenType == "hap")  // 只取单倍型的序列长度
        {
            if (gtVec.size() == 1)  // gt只有一个，代表为纯合的变异，两个单倍型添加一样的长度
            {
                if (gtVec[0] == 0)  // 如果基因型为0，跳过该位点
                {
                    // 返回 '0/0'
                    vector<int> qryLenVecTmp(seqLenVec[0], qrySeqVec.size());
                    return make_tuple(seqLenVec[0], qryLenVecTmp);
                }
                else
                {
                    qryLenVec.push_back(seqLenVec[gtVec[0]]);
                    qryLenVec.push_back(seqLenVec[gtVec[0]]);
                }
            }
            else
            {
                for (auto gtTmp : gtVec)
                {
                    qryLenVec.push_back(seqLenVec[gtTmp]);
                }
            }
        }
        else  // 取所有的长度
        {
            for (size_t i = 1; i < seqLenVec.size(); i++)
            {
                qryLenVec.push_back(seqLenVec[i]);
            }
        }

        return make_tuple(seqLenVec[0], qryLenVec);
    }


    /**
	 * 保存recall中failCall的结果
	 *
	 * @param chrStartLenInfoGtVecMap     真集中用剩下的vcf信息
     * @param outFileName                 输出文件名
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        const string & outFileName
    )
    {
        // 没有找到的变异
        string failCallTxt;
        // 输出文件流
        gzFile gzfp = gzopen(outFileName.c_str(), "wb");

        // 遍历字典，将没有找到的vcf进行保存
        for (auto it1 : chrStartLenInfoGtVecMap)  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        {
            for (auto it2 : it1.second)  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            {
                failCallTxt += get<2>(it2.second) + "\n";

                // 保存结果
                if (failCallTxt.size() > 10000000) // 每10Mb写入一次
                {
                    gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

                    // 清空字符串
                    failCallTxt.clear();
                    string().swap(failCallTxt);
                }
            }
        }
        // 保存结果
        gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

        // 释放内存，关闭文件
        gzclose(gzfp);

        return 0;
    }


    /**
	 * 保存recall中failCall的结果
	 *
	 * @param rocCallMap
     * @param rocRecallMap
     * @param all_length_list     所有变异长度
     * 
     * 
     * @return int             0
	**/
    void roc_calculate(
        const map<float, vector<int>> & rocCallMap, 
        const map<float, vector<int>> & rocRecallMap, 
        const vector<int> & all_length_list
    )
    {
        // 保存结果
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
        string outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        int allTruePositives = all_length_list.size(); // 所有正确的位点数量
        int truePositives = 0; // 真阳性
        int falsePositives = 0; // 假阳性
        int trueNegatives = 0; // 真阴性
        int falseNegatives = 0; // 假阴性

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            auto findIter1 = rocRecallMap.find(iter1->first);
            if (findIter1 == rocRecallMap.end())
            {
                falsePositives += iter1->second.size();
            }
            else
            {
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
        outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        allTruePositives = 0; // 所有正确的位点数量
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] == 0)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // 真阳性
        falsePositives = 0; // 假阳性
        trueNegatives = 0; // 真阴性
        falseNegatives = 0; // 假阴性

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all中的snp数量
            // 遍历该score下的长度vector
            for (size_t i = 0; i < iter1->second.size(); i++)
            {
                if (iter1->second[i] == 0)
                {
                    allPositiveTmp++;
                }
            }

            int truePositiveTmp = 0;
            auto findIter1 = rocRecallMap.find(iter1->first);
            if (findIter1 != rocRecallMap.end())
            {
                // 遍历该score下的长度vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] == 0)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // 如果个数是0，跳过该score
            if (allPositiveTmp == 0)
            {
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
        outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        allTruePositives = 0; // 所有正确的位点数量
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] > -50 && all_length_list[i] < 50 && all_length_list[i] != 0)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // 真阳性
        falsePositives = 0; // 假阳性
        trueNegatives = 0; // 真阴性
        falseNegatives = 0; // 假阴性

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all中的snp数量
            // 遍历该score下的长度vector
            for (size_t i = 0; i < iter1->second.size(); i++)
            {
                if (iter1->second[i] > -50 && iter1->second[i] < 50 && iter1->second[i] != 0)
                {
                    allPositiveTmp++;
                }
            }

            int truePositiveTmp = 0;
            auto findIter1 = rocRecallMap.find(iter1->first);
            if (findIter1 != rocRecallMap.end())
            {
                // 遍历该score下的长度vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // 如果个数是0，跳过该score
            if (allPositiveTmp == 0)
            {
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
        outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        allTruePositives = 0; // 所有正确的位点数量
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] <= -50)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // 真阳性
        falsePositives = 0; // 假阳性
        trueNegatives = 0; // 真阴性
        falseNegatives = 0; // 假阴性

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all中的snp数量
            // 遍历该score下的长度vector
            for (size_t i = 0; i < iter1->second.size(); i++)
            {
                if (iter1->second[i] <= -50)
                {
                    allPositiveTmp++;
                }
            }

            int truePositiveTmp = 0;
            auto findIter1 = rocRecallMap.find(iter1->first);
            if (findIter1 != rocRecallMap.end())
            {
                // 遍历该score下的长度vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] <= -50)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // 如果个数是0，跳过该score
            if (allPositiveTmp == 0)
            {
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
        outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        allTruePositives = 0; // 所有正确的位点数量
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] >= 50)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // 真阳性
        falsePositives = 0; // 假阳性
        trueNegatives = 0; // 真阴性
        falseNegatives = 0; // 假阴性

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all中的snp数量
            // 遍历该score下的长度vector
            for (size_t i = 0; i < iter1->second.size(); i++)
            {
                if (iter1->second[i] >= 50)
                {
                    allPositiveTmp++;
                }
            }

            int truePositiveTmp = 0;
            auto findIter1 = rocRecallMap.find(iter1->first);
            if (findIter1 != rocRecallMap.end())
            {
                // 遍历该score下的长度vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] >= 50)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // 如果个数是0，跳过该score
            if (allPositiveTmp == 0)
            {
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
} // namespace RECALL

#endif