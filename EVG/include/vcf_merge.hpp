#ifndef vcf_merge_hpp
#define vcf_merge_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <malloc.h>
#include "zlib.h"
#include <cmath>
#include "strip_split_join.hpp"
#include "sort_find.hpp"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"

using namespace std;

namespace MERGE
{ 
    struct baseVcfStruct
    {
        string headInfo; // vcf的注释行

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int>>> qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, 一个位点有多个等位基因的时候，用这个来保存qry的长度

        // 存储baseVcf的vcfInfo
        map<string, map<int, string>> BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

        // 存储二分查找法后的结果，保存软件名字和对应的分型结果
        map<string, map<int, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
        
        map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>>  最后两个为方差和标准差
        map<string, int> softwateSampleIdx;  // map<software, index>, 软件文件对应的sample的索引

        // 记录列数，用于Paragraph提取时，从该列之后找元素
        int colNum;
    };

    struct softwareVcfStruct
    {
        string software;  // 软件名字

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >，可能是1/2等情况，因此要存储每个等位基因的长度，如果长度为1代表为纯合变异
        map<string, vector<vector<int>>> gtVecMap;  // map<chromosome, vector<vector<genotype>>>
        map<string, vector<vector<float>>> depthVecMap;  // map<chromosome, vector<vector<depth>>>

        
        vector<float> depthVec;  // vector<depth>
    };

    /********************************************************************************************/
    /*
        检查文件是否存在
        fileName -> 需要检查的文件
    */
    bool check_file(
        string fileName
    );
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const int & sampleIdx
    );
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx = 0
    );
    int build_softwarefile_index(
        baseVcfStruct & mergeVcfStruct, 
        const string & vcfFileName, 
        const string &software, 
        const string & sample_name
    );
    baseVcfStruct build_basefile_index(
        const string & vcfFileName, 
        const string & sample_name
    );
    tuple<int, int> recall_push(
        baseVcfStruct & mergeVcfStruct, 
        const string chromosome, 
        const int trueRefStart, 
        const string software, 
        const vector<int> qryLenVec, 
        const vector<int> gtVec, 
        const int gt, 
        const float depth, 
        int j, 
        int k,
        int l,
        int & indexLeft, 
        int & indexRight
    );
    int vcf_merge(
        baseVcfStruct & mergeVcfStruct, 
        softwareVcfStruct & softvcfStructure
    );
    tuple<int, vector<int> > get_hap_len(
        const string & svType, 
        const string & refSeq, 
        const string & qrySeqs, 
        const vector<int> & gtVec, 
        const string & lenType
    );
    int sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    );
    string filter_rule(
        vector<vector<int>> & gtAllVec, 
        vector<string> & softwareVec, 
        vector<float> & depthNorVec, 
        const string & vcfInfo, 
        const int & softwareNum, 
        const string & mode, 
        baseVcfStruct & mergeVcfStruct
    );
    tuple<float, float, float> cal_var_sd(
        const vector<float> & data
    );
    int result_save(
        string & headInfo, 
        map<string, map<int, string > > & outMap, 
        const string & outputFileName
    );
    /********************************************************************************************/
    

   	/**
	 * Building index for different software.
	 *
	 * @param trueVcf              位点所有的gt的vector
	 * @param ParagraphVcf         Paragraph软件输出
	 * @param GraphTyper2Vcf       GraphTyper2软件输出
	 * @param BayesTyperVcf        BayesTyper软件输出
     * @param MAPVcf               VG-MAP软件输出
     * @param GiraffeVcf           VG-Giraffe软件输出
     * @param GraphAlignerVcf      GraphAligner软件输出
     * @param PanGenieVcf          PanGenie软件输出
     * @param mergeVcfStruct       合并后的struct
     * 
     * 
     * @return softwareNum
	**/

    int run_index_merge(
        string trueVcf, 
        string ParagraphVcf, 
        string GraphTyper2Vcf, 
        string BayesTyperVcf, 
        string MAPVcf, 
        string GiraffeVcf, 
        string GraphAlignerVcf, 
        string & PanGenieVcf, 
        baseVcfStruct & mergeVcfStruct, 
        string sample_name
    )
    {
        int softwareNum = 0;  // 记录软件的数量，用于对结果过滤

        if (trueVcf.length() > 0)
        {
            // 检查vcf文件是否排序
            check_vcf_sort(trueVcf);

            mergeVcfStruct = build_basefile_index(trueVcf, sample_name);
        }
        else
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Please enter the -v/--vcf parameter" << endl;
            exit(1);
        }

        vector<string> softwareVec = {ParagraphVcf, MAPVcf, GraphTyper2Vcf, BayesTyperVcf, GiraffeVcf, GraphAlignerVcf, PanGenieVcf};
        vector<string> softwareNameVec = {"Paragraph", "VG-MAP", "GraphTyper2", "BayesTyper", "VG-Giraffe", "GraphAligner", "PanGenie"};

        for (size_t i = 0; i < softwareVec.size(); i++)
        {
            if (softwareVec[i].length() > 0)
            {
                check_vcf_sort(softwareVec[i]);  // 检查vcf文件是否排序
                
                // 为软件输出特异的索引
                build_softwarefile_index(
                    mergeVcfStruct, 
                    softwareVec[i], 
                    softwareNameVec[i], 
                    sample_name
                );

                softwareNum++;  // 软件数量加1
            }
        }

        return softwareNum;
    }


    /*
        检查文件是否存在
        fileName -> 需要检查的文件
    */
    bool check_file(string fileName)
    {
        // 输入文件流
        gzFile gzfp = gzopen(fileName.c_str(), "rb");

        // 如果文件没有打开，则报错并退出程序。
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << fileName << "': No such file or directory." << endl;
            exit(1);
        }

        // 释放内存，关闭文件
        gzclose(gzfp);
        
        return true;
    }



    // paragraph DP -> 1
    // graphaligner, vg-map- vg-giraffe DP -> 6
    // bayestyper MAC -> -1,1.87332
    // graphtyper DP -> 6
    // pangenie KC -> 4
    /**
	 * 获取位点DP.
	 *
	 * @param informationsVec     vcfInfoVec
     * @param sampleIdx           sample对应的列
     * 
     * 
     * @return depthVec        vector<depth>
	**/
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const int & sampleIdx
    )
    {
        vector<float> depthVec;  // 位点覆盖度的vector

        vector<string> formatVec;  // FORMAT字段拆分
        int formatIndex = 8;  // FORMAT的索引
        formatVec = split(informationsVec[formatIndex], ":");

        int depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "DP"));  // 获取depth的索引位置

        if (depthIndex == formatVec.size())  // BayesTyper软件没有DP字段，所以用MAC字段来代替
        {
            depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "MAC"));
        }

        if (depthIndex == formatVec.size())  // PanGenie软件没有DP字段，所以用KC字段来代替
        {
            depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "KC"));
        }
        

        if (depthIndex == formatVec.size())  // 判断index是不是存在，不存在的话返回深度都为0。
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [DP,MAC,KC] information in FORMAT -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
            depthVec = {0, 0};
        }
        else  // 如果存在，则进行保存
        {
            for (auto it : split(split(informationsVec[sampleIdx], ":")[depthIndex], ","))  // 找到depth后，对其按，拆分并循环
            {
                if (it == "-1")  // BayesTyper的-1代表没找到，因此把-1改为0
                {
                    it = "0";
                }
                depthVec.push_back(stof(it));  // 添加到vector中
            }
        }

        return depthVec;
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
	 * 构件base文件索引，不对vcf进行过滤，用于创建基准index。文件为真实的vcf文件
	 *
	 * @param vcfFileName     base的vcf文件
	 * @param sample_name     样本名
     * 
     * 
     * @return outStruct
	**/
    baseVcfStruct build_basefile_index(
        const string & vcfFileName, 
        const string & sample_name
    )
    {
        baseVcfStruct outStruct;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: '" << vcfFileName << "'" << endl;

        vector<int> all_length_list;  // 循环将vcf信息加到outStructure中

        gzFile gzfp = gzopen(vcfFileName.c_str(), "rb");  // 输入文件流

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "'" << vcfFileName << "': No such file or directory." << endl;
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

                    string::size_type findIdx1 = informations.find("#");  // 注释行索引

                    if (findIdx1 != string::npos)  // 注释行
                    {
                        string::size_type findIdx2 = informations.find("#CHROM");  // #CHROM注释行索引
                        if (findIdx2 != string::npos)
                        {
                            vector<string> informationsVec = split(informations, "\t");  // vcfInfo拆分
                            vector<string> informationsVecTmp{&informationsVec[0], &informationsVec[0]+9};  // 将FORMAT后边的信息全部去掉，然后改为 sample_name
                            outStruct.headInfo += join(informationsVecTmp, "\t") + "\t" + sample_name + "\n";

                            // 记录baseVcf的列数
                            outStruct.colNum = informationsVec.size();
                        }
                        else
                        {
                            outStruct.headInfo += informations + "\n";
                        }
                    }
                    else
                    {
                        
                        int qryLen;
                        
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo拆分
                        string chromosome = informationsVec[0];  // 染色体信息
                        int refStart = stoi(informationsVec[1]);  // 变异的起始位置
                        string svType = informationsVec[2];  // 判断BayesTyper的结果为duplication
                        string refSeq = informationsVec[3];  // 变异的ref序列
                        string qrySeqs = informationsVec[4];  // 变异的qry序列

                        outStruct.startMap[chromosome].push_back(refStart);  // 保存变异起始位置

                        vector<string> informationsVecTmp{&informationsVec[0], &informationsVec[0]+8};  // 将FORMAT改成GT，后边的信息全部去掉
                        outStruct.BaseInfoMap[chromosome][refStart] = join(informationsVecTmp, "\t") + "\tGT:SO:DP:ADP:NDP";  // 保存变异informationsVecTmp信息 GT:software:depth:averageDepth:normalDepth

                        // 获取长度信息
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            {0, 0}, 
                            "all"
                        );

                        // 保存变异的ref长度和qry长度。qry为vector，因为有多个等位基因。
                        outStruct.refLenMap[chromosome].push_back(refLen);
                        outStruct.qryLenVecMap[chromosome].push_back(qryLenVec);

                        //清空vector，释放内存
                        vector<int>().swap(qryLenVec);
                    }
                    // 清空字符串
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // 释放内存，关闭文件
        gzclose(gzfp);

        return outStruct;
    }


    /**
	 * 对软件的结果进行过滤并构建索引，文件为软件的输出文件。
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index构建的base索引
	 * @param vcfFileName        软件输出的vcf文件
     * @param software           软件名
     * @param sample_name        样本名，用于提取基因型信息
     * 
     * 
     * @return 0
	**/
    int build_softwarefile_index(
        baseVcfStruct & mergeVcfStruct, 
        const string & vcfFileName, 
        const string & software, 
        const string & sample_name
    )
    {
        softwareVcfStruct softvcfStructure;  // 输出struct

        softvcfStructure.software = software;  // 软件名

        int sampleIdx = 0;  // 记录sample在vcf文件的第几列

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: " << software << endl;  // print log

        vector<int> all_length_list;  // 循环将vcf信息加到outStructure中

        gzFile gzfp = gzopen(vcfFileName.c_str(), "rb");  // 输入文件流

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << vcfFileName << "': No such file or directory." << endl;
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

                    string::size_type findidx1 = informations.find("#");  // 注释行索引
                    if (findidx1 != string::npos)  // 注释行
                    {
                        string::size_type findIdx2 = informations.find("#CHROM");  // #CHROM注释行索引
                        if (findIdx2 != string::npos)  // 如果是的话找sample_name索引
                        {
                            vector<string> informationsVec = split(informations, "\t");  // vcfInfo拆分

                            // 获取sample对应的列数
                            if (software == "Paragraph")  // 因为paragraph会保留原始的列，所以从原始列之后找sample  mergeVcfStruct.colNum
                            {
                                sampleIdx = distance(informationsVec.begin(), find(informationsVec.begin()+mergeVcfStruct.colNum, informationsVec.end(), sample_name));  // 获取sample的索引位置
                            }
                            else  // 其它软件都从第10列后找 FORMAT在第9列 
                            {
                                sampleIdx = distance(informationsVec.begin(), find(informationsVec.begin()+9, informationsVec.end(), sample_name));  // 获取sample的索引位置
                            }
                            
                            // 如果没找到，报错并退出代码
                            if (sampleIdx == informationsVec.size())
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " 
                                     << "Error: '" << sample_name << "' is not in the output file of " 
                                     << software 
                                     << " -> " 
                                     << informations 
                                     << endl;
                                exit(1);
                            }

                            mergeVcfStruct.softwateSampleIdx[software] = sampleIdx; // 软件对应的sample索引添加
                        }
                    }
                    else  // 非注释行
                    {
                        // 判断sample_name找到没有
                        if (sampleIdx == 0)
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " 
                                << "Error: '" << sample_name << "' is not in the output file of " 
                                << software 
                                << endl;
                            exit(1);
                        }

                        vector<string> informationsVec = split(informations, "\t");
                        string chromosome = informationsVec[0];
                        int refStart = stoi(informationsVec[1]);
                        string svType = informationsVec[2];  // 判断BayesTyper的结果为duplication
                        string refSeq = informationsVec[3];
                        string qrySeqs = informationsVec[4];
                        string filter = informationsVec[6];
                        
                        // 获取位点的覆盖度信息
                        vector<float> depthVec = get_depth(
                            informationsVec, 
                            sampleIdx
                        );
                        float depSum = accumulate(depthVec.begin(), depthVec.end(), 0);  // 计算所有等位基因的深度和
                        softvcfStructure.depthVec.push_back(depSum);  // 软件对应的vcf深度列表添加
                        
                        vector<int> gtVec = get_gt(
                            informationsVec, 
                            sampleIdx
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

                        // 根据FILTER字段进行过滤
                        if (filter != "PASS") //如果不是(PASS)则直接跳过，不构建索引
                        {
                            // 清空字符串
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                        
                        // 检查索引是否越界
                        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
                        if (maxGtNum > split(qrySeqs, ",").size())  // 先检查是否数组是否越界
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " 
                                << "Error: number of genotyping and query sequences do not match -> " 
                                << informations << endl;
                            exit(1);
                        }

                        // 获取单倍型的长度信息
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // 保存变异的起始和分型信息
                        softvcfStructure.startMap[chromosome].push_back(refStart);
                        softvcfStructure.gtVecMap[chromosome].push_back(gtVec);
                        softvcfStructure.depthVecMap[chromosome].push_back(depthVec);  // 软件对应的vcf深度（每个等位基因的）

                        // 保存变异的refLen和qryLenVec
                        softvcfStructure.refLenMap[chromosome].push_back(refLen);
                        softvcfStructure.qryLenVecMap[chromosome].push_back(qryLenVec);
                    }
                    // 清空字符串
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // 释放内存，关闭文件
        gzclose(gzfp);

        // 计算软件的平均覆盖度，方差和标准差
        float aveDepth;
        float variance;
        float sd;
        tie(aveDepth, variance, sd) = cal_var_sd(softvcfStructure.depthVec);
        get<0>(mergeVcfStruct.depthMap[software]) = aveDepth;
        get<1>(mergeVcfStruct.depthMap[software]) = variance;
        get<2>(mergeVcfStruct.depthMap[software]) = sd;

        // 将软件的结果加到图中
        vcf_merge(mergeVcfStruct, softvcfStructure);

        return 0;
    }


    /**
	 * vcf_merge中recall后添加到总图中函数
	 *
     * @date 2023/3/7
	 * @param mergeVcfStruct     build_basefile_index构建的base索引
	 * @param chromosome         染色体号
     * @param trueRefStart       真实的ref起始位置
     * @param software           软件
     * @param qryLenVec          软件分型的单倍型对应长度
     * @param gtVec              软件的分型gt列表
     * @param gt                 当前循环的gt
     * @param depth              位点对应的深度
     * @param j                  软件多个单倍型循环的索引
     * @param k                  二分查找法两个坐标循环的索引
     * @param l                  真实位点单倍型的循环索引
     * @param indexLeft          二分查找法的左索引
     * @param indexRight         二分查找法的右索引
     * 
     * 
     * @return tuple<int, int>   tuple<state, j>  0 -> 正确添加；-1 -> 需要进行下一个循环
	**/
    tuple<int, int> recall_push(
        baseVcfStruct & mergeVcfStruct, 
        const string chromosome, 
        const int trueRefStart, 
        const string software, 
        const vector<int> qryLenVec, 
        const vector<int> gtVec, 
        const int gt, 
        const float depth, 
        int j, 
        int k,
        int l,
        int & indexLeft, 
        int & indexRight
    )
    {
        // 先初始化recall哈希表
        if (mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome].find(trueRefStart) == mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome].end())  // 如果起始位置第一次见，则初始化
        {
            mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart];
        }
        if (mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart].find(software) == mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart].end())  // 如果染色体第一次见，则初始化
        {
            vector<float> depthVecTmp(qryLenVec.size());
            vector<int> gtVecTmp(qryLenVec.size());
            mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software] = make_tuple(depthVecTmp, gtVecTmp);
        }

        // gt和l中有一个时0，但另一个不是，则下一个循环，防止SNP判断时候出错
        if ((gt != 0 && l == 0) || (gt == 0 && l != 0))
        {
            return make_tuple(-1, j);
        }
        else
        {
            // 每个单倍型单独添加
            get<0>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;  // 对应单倍型位置深度
            get<1>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;  // 对应单倍型位置基因型

            // 更改二分查找的坐标，下一个等位基因就直接添加到该refStart上
            indexLeft = k; 
            indexRight = k;

            if (j == 0 && gt == gtVec[1])  // 如果两个单倍型是一样的，直接赋值，然后跳过下一个循环
            {
                j++;
                get<0>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;  // 对应单倍型位置深度
                get<1>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;  // 对应单倍型位置基因型
            }
        }

        return make_tuple(0, j);
    }


    /**
	 * 对软件的结果进行过滤并构建索引，文件为软件的输出文件。
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index构建的base索引
	 * @param softvcfStructure   build_softwarefile_index构建的软件vcf索引
     * 
     * 
     * @return 0
	**/
    int vcf_merge(
        baseVcfStruct & mergeVcfStruct, 
        softwareVcfStruct & softvcfStructure
    )
    {
        vector<string> chromosomeVec;
        string software = softvcfStructure.software;

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Merging: " << software << ".\n";  // print log

        for (auto iter = mergeVcfStruct.startMap.begin(); iter != mergeVcfStruct.startMap.end(); iter++)
        {
            chromosomeVec.push_back(iter->first);
        }

        for (auto iter = softvcfStructure.startMap.begin(); iter != softvcfStructure.startMap.end(); iter++)  // 循环软件的vcf（染色体）
        {
            string chromosome = iter->first;
            vector<int> startVec = softvcfStructure.startMap[chromosome];
            vector<int> refLenVec = softvcfStructure.refLenMap[chromosome];
            vector<vector<int> > qryLenVecVec = softvcfStructure.qryLenVecMap[chromosome];
            vector<vector<int> > gtVecVec = softvcfStructure.gtVecMap[chromosome];
            vector<vector<float> > depthVecVec = softvcfStructure.depthVecMap[chromosome];

            // 先初始化recall哈希表
            if (mergeVcfStruct.recallSoftwareGtDepVecMap.find(chromosome) == mergeVcfStruct.recallSoftwareGtDepVecMap.end())  // 如果染色体第一次见，则初始化
            {
                mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome];
            }

            if (find(chromosomeVec.begin(), chromosomeVec.end(), chromosome) == chromosomeVec.end())  // 先看mergeVcfStructure中有没有该染色体，没有的话报错并退出代码
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: chromosome is not present in the base vcf file -> " << chromosome << endl;
                exit(1);
            }
            else  // 如果存在的话，则用二分查找法找最近的变异，判断它们是不是同一个，是的话push，不是的话添加一个新的
            {
                vector<int> mergeRefStartVec = mergeVcfStruct.startMap[chromosome];
                vector<int> mergeRefLenVec = mergeVcfStruct.refLenMap[chromosome];
                vector<vector<int> > mergeQryLenVecVec = mergeVcfStruct.qryLenVecMap[chromosome];  // vector<vector<qryLen>>

                for (int i = 0; i < startVec.size(); i++)  // 循环软件的vcf（位置）
                {
                    int refStart = startVec[i];
                    int refLen = refLenVec[i];
                    vector<int> qryLenVec = qryLenVecVec[i];  // 位点分型对应的等位基因序列长度
                    vector<int> gtVec = gtVecVec[i];  // 位点分型对应的分型
                    vector<float> depthVec = depthVecVec[i];  // 位点分型对应的等位基因深度

                    // 记录二分查找法的索引，多个等位基因用同一个refStart
                    int indexLeft = -1;
                    int indexRight = -1;

                    // 对多个等位基因遍历
                    for (size_t j = 0; j < qryLenVec.size(); j++)
                    {
                        int qryLen = qryLenVec[j];
                        int gt = gtVec[j];
                        float depth = depthVec[j];
                    
                        // 二分查找法找200/10bp内的变异索引。
                        if (indexLeft == -1 && indexRight == -1)  // 新的等位基因再查找
                        {
                            if (refLen <= 49 && qryLen <= 49)
                            {
                                indexLeft = search_Binary_left(mergeRefStartVec, (refStart-10));
                                indexRight = search_Binary_right(mergeRefStartVec, (refStart+10));
                            }
                            else
                            {
                                indexLeft = search_Binary_left(mergeRefStartVec, (refStart-200));
                                indexRight = search_Binary_right(mergeRefStartVec, (refStart+200));
                            }

                            if (indexLeft < 0 || indexRight >= mergeRefStartVec.size())
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: out of index, please check the data or code.\n";
                                exit(1);
                            }
                        }

                        for (int k = indexLeft; k <= indexRight; k++)
                        {
                            // true变异的信息
                            vector<int> mergeQryLenVec = mergeQryLenVecVec[k];  // 多个等位基因的长度
                            int trueRefLen = mergeRefLenVec[k];
                            int trueRefStart = mergeRefStartVec[k];

                            // 构造ref和qry共同的长度索引
                            vector<int> trueSeqLenVec;
                            trueSeqLenVec.push_back(trueRefLen);  // 先添加参考基因组的长度
                            trueSeqLenVec.insert(trueSeqLenVec.end(), mergeQryLenVec.begin(), mergeQryLenVec.end());  // 再添加

                            // 如果长度为0，则报错，并退出代码。
                            if (mergeQryLenVec.size() == 0)
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " 
                                    << "Error: mergeQryLenVec.size() == 0 -> chromosome:" << chromosome 
                                    << " refStart:" << mergeRefStartVec[k] << endl;
                                exit(1);
                            }

                            for (size_t l = 0; l < trueSeqLenVec.size(); l++)  // 一个位点有多个等位基因的时候，trueSeqLenVec存储了各个等位基因的长度，因此遍历它看该位点的变异是否和merge里边的一致，包含0
                            {
                                int trueQryLen = trueSeqLenVec[l];

                                // 缺失
                                if (refLen >= 50 && qryLen < 50)
                                {
                                    if ((abs(refStart-trueRefStart)<=200) &&
                                        (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                        ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
                                    {
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
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
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }

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
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
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
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
                                        goto stop;
                                    }
                                }
                                // del
                                else if (refLen >= 3 && qryLen <= 2)
                                {
                                    if ((abs(refStart-trueRefStart)<=1) &&
                                    (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                    ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
                                    {
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
                                        goto stop;
                                    }
                                }
                                // ins
                                else if (refLen <= 2 && qryLen >= 3)
                                {
                                    if ((abs(refStart-trueRefStart)<=1) &&
                                    (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                    ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
                                    {
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
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
                                        // 添加到总图中
                                        int recall_state;
                                        tie(recall_state, j) = recall_push(
                                            mergeVcfStruct, 
                                            chromosome, 
                                            trueRefStart, 
                                            software, 
                                            qryLenVec, 
                                            gtVec, 
                                            gt, 
                                            depth, 
                                            j, 
                                            k,
                                            l,
                                            indexLeft, 
                                            indexRight
                                        );

                                        // 没有添加则下一个循环
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
                                        goto stop;
                                    }
                                }
                            }
                        }
                        stop:; // 如果找到了，则退出嵌套循环。结束该单倍型的循环，继续下一个单倍型
                    }
                }
            }
        }
        return 0;
    }


    /**
	 * 获取变异的长度信息
	 *
     * @param vcfInfo           vcfInfo
     * @param gtVec             基因型列表
     * 
     * 
     * @return int              svLength
	**/
    int sv_length_select(
        const string & vcfInfo, 
        const vector<int> & gtVec
    )
    {
        int svLength;

        vector<string> vcfInfoVec = split(vcfInfo, "\t");

        // 获取ref和单倍型的长度
        int refLen;
        vector<int> qryLenVec;
        tie(refLen, qryLenVec) = get_hap_len(
            vcfInfoVec[2], 
            vcfInfoVec[3], 
            vcfInfoVec[4], 
            gtVec, 
            "hap"
        );

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
	 * Filter result according to the depth and number of variants.
	 *
	 * @param gtAllVec        位点所有的gt的vector
	 * @param softwareVec     位点所有的software的vector
	 * @param depthNorVec     位点所有标准化后的depth的vector
     * @param softwareNum     用户输入的软件数量
     * @param vcfInfo         位点的变异信息
     * @param mode            EVG运行的模式
     * @param mergeVcfStruct  总图
     * 
     * @return software
	**/
    string filter_rule(
        vector<vector<int>> & gtAllVec, 
        vector<string> & softwareVec, 
        vector<float> & depthNorVec, 
        const string & vcfInfo, 
        const int & softwareNum, 
        const string & mode, 
        baseVcfStruct & mergeVcfStruct
    )
    {
        // 结果过滤
        vector<string> bestSoftwareVec{};  // 出现频率最高的软件vector
        string selectSoftware;  // 根据sv的长度进行筛选后的软件vector
        int svLength{};
        int gtFrequency = 0;

        // 标准化后的深度最接近0的索引
        float selectDepth = 10000.0;
        int selectDepthIdx = -1;
        for (int i = 0; i < depthNorVec.size(); i++)
        {
            if (abs(depthNorVec[i]) < selectDepth)
            {
                selectDepth = abs(depthNorVec[i]);  // 更新深度
                selectDepthIdx = i;  // 更新索引
            }
        }

        // 计算每一种基因型出现的频率
        // 先找频率最大的gt
        for (size_t i = 0; i < gtAllVec.size(); i++)
        {
            svLength = sv_length_select(
                vcfInfo, 
                gtAllVec[i]
            );

            if (count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]) > gtFrequency)
            {
                gtFrequency = count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]);
            }
        }

        // 把频率最大的gt对应软件找出来
        for (size_t i = 0; i < gtAllVec.size(); i++)
        {
            if (count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]) == gtFrequency)
            {
                bestSoftwareVec.push_back(softwareVec[i]);
            }
        }

        // 过滤规则
        if (-49 <= svLength && svLength <= 49)  // snp+indel
        {
            if (mode == "specific")  // 为株系特异的vcf文件
            {
                if (gtFrequency < 2) // 支持的软件数量小于2，选取深度最接近于0的基因型
                {
                    selectSoftware = softwareVec[selectDepthIdx];  // 选取覆盖度接近于0的软件
                }
                else // 结果有两个以上软件支持的话，则选择该基因型。
                {
                    selectSoftware = bestSoftwareVec[0]; // 选择出现频率最高的分型结果
                }
            }
            else  // 为构图用的vcf
            {
                if (gtFrequency < 2) // 结果必须有两个以上的软件支持，否则跳过该位点
                {
                    selectSoftware.clear();
                }
                else // 结果有两个以上软件支持的话，则选择该基因型。
                {
                    selectSoftware = bestSoftwareVec[0]; // 选择出现频率最高的分型结果
                }
            }
        }
        else  // Insertion and Deletion
        {
            vector<string> softwareSortVec = {"BayesTyper", "Paragraph", "VG-MAP", "VG-Giraffe", "GraphTyper2", "GraphAligner", "PanGenie"};
            for (auto it : softwareSortVec)
            {
                // 大于50bp的变异中如果有it支持，则直接选择it的结果
                if (find(softwareVec.begin(), softwareVec.end(), it) != softwareVec.end())
                {
                    selectSoftware = it;
                    break;
                }
                else // 否则选取覆盖度接近0的软件
                {
                    selectSoftware = softwareVec[selectDepthIdx];  // 选取覆盖度接近于0的软件
                }
            }
        }

        return selectSoftware;
    }



    /*
        vcf过滤
        mergeVcfStruct -> 软件合并后的struct
    */

   /**
	 * Filter and merge result according to the depth and number of variants.
	 *
	 * @param mergeVcfStruct     软件合并后的struct
	 * @param softwareNum        用户提交的软件数量
	 * @param mode               EVG运行的模式
     * 
     * @return outMap            outMap<chromosome, map<refStart, vcfInfo> >
	**/
   map<string, map<int, string > > vcf_merge_filter(
        baseVcfStruct & mergeVcfStruct, 
        const int & softwareNum, 
        const string & mode
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Filtering.\n";

        // 保存结果
        map<string, map<int, string > > outMap;  // outMap<chromosome, map<refStart, vcfInfo> >

        // 遍历recall哈希表
        for (auto it1 : mergeVcfStruct.recallSoftwareGtDepVecMap)  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
        {
            string chromosome = it1.first;
            
            for (auto it2 : it1.second) // map<start, map<software, tuple<vector<depth>, vector<gt> > > >
            {
                int refStart = it2.first;

                // 该位点的vcfInfo
                string informations = mergeVcfStruct.BaseInfoMap[chromosome][refStart];

                // 构造filter_rule的参数
                vector<vector<int>> gtAllVec;
                vector<string> softwareVec;
                vector<float> depthSumVec;
                vector<float> depthNorVec;
                vector<vector<float> > depthVecVec;

                for (auto it3 : it2.second)  // map<software, tuple<vector<depth>, vector<gt> > >
                {
                    string softwareTmp = it3.first;
                    vector<float> depthVecTmp = get<0>(it3.second);
                    vector<int> gtVecTmp = get<1>(it3.second);

                    depthVecVec.push_back(depthVecTmp);

                    // 获取软件平均的覆盖度，先检查map中有没有，如果没有则报错
                    float meanDepth;
                    float sd;
                    if (mergeVcfStruct.depthMap.find(softwareTmp) != mergeVcfStruct.depthMap.end())
                    {
                        meanDepth = get<0>(mergeVcfStruct.depthMap[softwareTmp]);
                        sd = get<2>(mergeVcfStruct.depthMap[softwareTmp]);
                    }
                    else
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: no coverage information -> "<< softwareTmp << endl;
                        exit(1);
                    }
                    float depthSumTmp = accumulate(depthVecTmp.begin(), depthVecTmp.end(), 0);  // 计算所有等位基因的深度和
                    depthSumVec.push_back(depthSumTmp);
                    
                    // 对覆盖度标准化（z-score标准化）
                    float depthNor = (depthSumTmp - meanDepth) / (float)sd;
                    
                    gtAllVec.push_back(gtVecTmp);
                    softwareVec.push_back(softwareTmp);
                    depthNorVec.push_back(depthNor);
                }

                // 过滤，选择最有可能的单倍型
                string selectSoftware = filter_rule(
                    gtAllVec, 
                    softwareVec, 
                    depthNorVec, 
                    informations, 
                    softwareNum, 
                    mode, 
                    mergeVcfStruct
                );

                if (selectSoftware.length() == 0) // 如果返回的是空的字符串，则表明位点不可信，跳过该位点
                {
                    continue;
                }

                // 选择的软件索引、基因型和深度信息
                int softwareIdx = distance(softwareVec.begin(), find(softwareVec.begin(), softwareVec.end(), selectSoftware));  // 获取software的索引位置
                vector<int> selectGtVec = gtAllVec[softwareIdx];
                float selectDepthSumVec = depthSumVec[softwareIdx];
                float selectDepthNor = depthNorVec[softwareIdx];
                vector<float> selectDepthVec = depthVecVec[softwareIdx];

                // // 保存结果
                string selectGt = join(selectGtVec, "/");
                selectGt += ":" + selectSoftware 
                        + ":" + to_string(selectDepthNor).substr(0, 5) 
                        + ":" + to_string(get<0>(mergeVcfStruct.depthMap[selectSoftware])).substr(0, 5) 
                        + ":" + to_string(selectDepthSumVec).substr(0, 5);
                outMap[chromosome][refStart] = informations + "\t" + selectGt;
            }
        }

        return outMap;
    }



    /**
	 * 计算方差和标准差
	 *
	 * @param data     含有位点深度的列表
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
	**/
    tuple<float, float, float> cal_var_sd(const vector<float>& data)
    {
        float sum = std::accumulate(std::begin(data), std::end(data), 0.0);
        float mean = sum / data.size();

        float variance = 0.0;
        std::for_each(std::begin(data), std::end(data), [&](const float d) {
            variance += pow(d-mean, 2);
        });
        variance /= data.size();

        float std_deviation = sqrt(variance);

        return make_tuple(mean, variance, std_deviation);
    }


    /*
        保存结果
        mergeVcfStruct -> 软件合并后的struct
        outMap -> vcf_merge_filter输出结果   outMap[chromosome][refStart][refLen][qryLen]
        prefix -> 输出文件名前缀
    */

    /**
	 * 保存结果
	 *
     * @param headInfo           vcf注释行
	 * @param outMap             vcf_merge_filter输出结果   outMap[chromosome][refStart]
	 * @param outputFileName     输出文件名
     * 
     * @return tuple   make_tuple(mean, variance, std_deviation)
	**/
    int result_save(
        string & headInfo, 
        map<string, map<int, string > > & outMap, 
        const string & outputFileName
    )
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Save result.\n";

        // 保存需要输出的string
        string outTxt = headInfo;

        for (auto it1 : outMap)  // outMap<chromosome, map<refStart, vcfInfo> >
        {
            for (auto it2: it1.second)  // map<refStart, vcfInfo>
            {
                string vcfInfo = it2.second;
                outTxt += vcfInfo + "\n";
            }
        }

        if (outputFileName.size() == 0)  // 没有指定输出文件，打印到屏幕
        {
            cout << outTxt << endl;
        }
        else if (outputFileName.find(".gz") != string::npos)  // 输出为压缩文件
        {
            // 输出文件流
            gzFile gzfp = gzopen(outputFileName.c_str(), "wb");

            // 打开文件
            if(!gzfp)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'"
                    << outputFileName 
                    << "': No such file or directory." 
                    << endl;
                exit(1);
            }
            else
            {
                gzwrite(gzfp, outTxt.c_str(), outTxt.length());
            }
            // 释放内存，关闭文件
            gzclose(gzfp);
        }
        else
        {
            // 输出文件流
            ofstream outputFile;
            outputFile.open(outputFileName, ios::out);

            // 打开文件
            if(!outputFile)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'"
                    << outputFileName 
                    << "': No such file or directory." 
                    << endl;
                exit(1);
            }
            else
            {
                outputFile << outTxt << endl;
            }

            // 关闭文件
            outputFile.close();
        }

        return 0;
    }
}  // namespace MERGE

#endif