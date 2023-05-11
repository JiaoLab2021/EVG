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
        string headInfo; // vcf��ע����

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int>>> qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >, һ��λ���ж����λ�����ʱ�������������qry�ĳ���

        // �洢baseVcf��vcfInfo
        map<string, map<int, string>> BaseInfoMap;  // map<chromosome, map<start, vcfInfo>>

        // �洢���ֲ��ҷ���Ľ��������������ֺͶ�Ӧ�ķ��ͽ��
        map<string, map<int, map<string, tuple<vector<float>, vector<int> > > > > recallSoftwareGtDepVecMap;  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
        
        map<string, tuple<float, float, float>> depthMap;  // map<software, pair<aveDepth, variance, sd>>  �������Ϊ����ͱ�׼��
        map<string, int> softwateSampleIdx;  // map<software, index>, ����ļ���Ӧ��sample������

        // ��¼����������Paragraph��ȡʱ���Ӹ���֮����Ԫ��
        int colNum;
    };

    struct softwareVcfStruct
    {
        string software;  // �������

        map<string, vector<int>> startMap;  // map<chromosome, vector<start>>
        map<string, vector<int>> refLenMap;  // map<chromosome, vector<refLen>>
        map<string, vector<vector<int> > > qryLenVecMap;  // map<chromosome, vector<vector<qryLen> > >��������1/2����������Ҫ�洢ÿ����λ����ĳ��ȣ��������Ϊ1����Ϊ���ϱ���
        map<string, vector<vector<int>>> gtVecMap;  // map<chromosome, vector<vector<genotype>>>
        map<string, vector<vector<float>>> depthVecMap;  // map<chromosome, vector<vector<depth>>>

        
        vector<float> depthVec;  // vector<depth>
    };

    /********************************************************************************************/
    /*
        ����ļ��Ƿ����
        fileName -> ��Ҫ�����ļ�
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
	 * @param trueVcf              λ�����е�gt��vector
	 * @param ParagraphVcf         Paragraph������
	 * @param GraphTyper2Vcf       GraphTyper2������
	 * @param BayesTyperVcf        BayesTyper������
     * @param MAPVcf               VG-MAP������
     * @param GiraffeVcf           VG-Giraffe������
     * @param GraphAlignerVcf      GraphAligner������
     * @param PanGenieVcf          PanGenie������
     * @param mergeVcfStruct       �ϲ����struct
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
        int softwareNum = 0;  // ��¼��������������ڶԽ������

        if (trueVcf.length() > 0)
        {
            // ���vcf�ļ��Ƿ�����
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
                check_vcf_sort(softwareVec[i]);  // ���vcf�ļ��Ƿ�����
                
                // Ϊ���������������
                build_softwarefile_index(
                    mergeVcfStruct, 
                    softwareVec[i], 
                    softwareNameVec[i], 
                    sample_name
                );

                softwareNum++;  // ���������1
            }
        }

        return softwareNum;
    }


    /*
        ����ļ��Ƿ����
        fileName -> ��Ҫ�����ļ�
    */
    bool check_file(string fileName)
    {
        // �����ļ���
        gzFile gzfp = gzopen(fileName.c_str(), "rb");

        // ����ļ�û�д򿪣��򱨴��˳�����
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << fileName << "': No such file or directory." << endl;
            exit(1);
        }

        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);
        
        return true;
    }



    // paragraph DP -> 1
    // graphaligner, vg-map- vg-giraffe DP -> 6
    // bayestyper MAC -> -1,1.87332
    // graphtyper DP -> 6
    // pangenie KC -> 4
    /**
	 * ��ȡλ��DP.
	 *
	 * @param informationsVec     vcfInfoVec
     * @param sampleIdx           sample��Ӧ����
     * 
     * 
     * @return depthVec        vector<depth>
	**/
    vector<float> get_depth(
        const vector<string> & informationsVec, 
        const int & sampleIdx
    )
    {
        vector<float> depthVec;  // λ�㸲�Ƕȵ�vector

        vector<string> formatVec;  // FORMAT�ֶβ��
        int formatIndex = 8;  // FORMAT������
        formatVec = split(informationsVec[formatIndex], ":");

        int depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "DP"));  // ��ȡdepth������λ��

        if (depthIndex == formatVec.size())  // BayesTyper���û��DP�ֶΣ�������MAC�ֶ�������
        {
            depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "MAC"));
        }

        if (depthIndex == formatVec.size())  // PanGenie���û��DP�ֶΣ�������KC�ֶ�������
        {
            depthIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "KC"));
        }
        

        if (depthIndex == formatVec.size())  // �ж�index�ǲ��Ǵ��ڣ������ڵĻ�������ȶ�Ϊ0��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [DP,MAC,KC] information in FORMAT -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
            depthVec = {0, 0};
        }
        else  // ������ڣ�����б���
        {
            for (auto it : split(split(informationsVec[sampleIdx], ":")[depthIndex], ","))  // �ҵ�depth�󣬶��䰴����ֲ�ѭ��
            {
                if (it == "-1")  // BayesTyper��-1����û�ҵ�����˰�-1��Ϊ0
                {
                    it = "0";
                }
                depthVec.push_back(stof(it));  // ��ӵ�vector��
            }
        }

        return depthVec;
    }


    /**
	 * ��ȡλ��������б�.
	 *
	 * @param informationsVec  vcfInfoList
     * @param sampleIdx        sample�����͵�����,Ĭ��ֵ0�������һ��
     * 
     * 
     * @return gtVec           vector <int>
	**/
    vector<int> get_gt(
        const vector<string> & informationsVec, 
        int sampleIdx
    )
    {
        vector <int> gtVec;  // λ����͵�vector

        vector<string> formatVec;  // FORMAT�ַ����
        int formatIndex = 8; // FORMAT������
        formatVec = split(informationsVec[formatIndex], ":");

        int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // ��ȡGT������λ��

        if (gtIndex == formatVec.size())  // �ж�index�Ƿ���ڣ������ڵĻ����ػ����Ͷ�Ϊ0��
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << informationsVec[0] << ":" << informationsVec[1] << endl;
            gtVec = {0, 0};
        }
        else  // ������ڣ�����б���
        {
            string gt;  // �洢�������ֶ�

            if (sampleIdx == 0) // û��ָ�������������һ��
            {
                gt = split(informationsVec[informationsVec.size()-1], ":")[gtIndex];  // gt�ֶ�
            }
            else
            {
                gt = split(informationsVec[sampleIdx], ":")[gtIndex];  // gt�ֶ�
            }
            
            string splitStr;  // gt�еķָ���
            if (gt.find("/") != string::npos)  // �жϡ�/���ָ���
            {
                splitStr = "/";
            }
            else if (gt.find("|") != string::npos)  // �жϡ�|��Ϊ�ָ���
            {
                splitStr = "|";
            }
            else  // ��֪����ʱ��Ϊ���ؿ�ֵ
            {
                gtVec = {0, 0};
                return gtVec;
            }
            
            for (auto it : split(gt, splitStr))  // �ҵ�gt�󣬶��䰴splitStr��ֲ�ѭ��
            {
                if (it == ".")  // ���Ϊ'.'��������λ��
                {
                    gtVec = {0, 0};
                    return gtVec;
                }
                gtVec.push_back(stoi(it));  // ��ӵ�vector��
            }
        }

        return gtVec;
    }


    /**
	 * ����base�ļ�����������vcf���й��ˣ����ڴ�����׼index���ļ�Ϊ��ʵ��vcf�ļ�
	 *
	 * @param vcfFileName     base��vcf�ļ�
	 * @param sample_name     ������
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

        vector<int> all_length_list;  // ѭ����vcf��Ϣ�ӵ�outStructure��

        gzFile gzfp = gzopen(vcfFileName.c_str(), "rb");  // �����ļ���

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "'" << vcfFileName << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // һ��ֻ��1024�ֽڵ�����
            while(gzgets(gzfp, line, 1024))
            {
                informations += line;
                if (informations.find("\n") != string::npos) // һ�н���
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // ȥ�����з�

                    string::size_type findIdx1 = informations.find("#");  // ע��������

                    if (findIdx1 != string::npos)  // ע����
                    {
                        string::size_type findIdx2 = informations.find("#CHROM");  // #CHROMע��������
                        if (findIdx2 != string::npos)
                        {
                            vector<string> informationsVec = split(informations, "\t");  // vcfInfo���
                            vector<string> informationsVecTmp{&informationsVec[0], &informationsVec[0]+9};  // ��FORMAT��ߵ���Ϣȫ��ȥ����Ȼ���Ϊ sample_name
                            outStruct.headInfo += join(informationsVecTmp, "\t") + "\t" + sample_name + "\n";

                            // ��¼baseVcf������
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
                        
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo���
                        string chromosome = informationsVec[0];  // Ⱦɫ����Ϣ
                        int refStart = stoi(informationsVec[1]);  // �������ʼλ��
                        string svType = informationsVec[2];  // �ж�BayesTyper�Ľ��Ϊduplication
                        string refSeq = informationsVec[3];  // �����ref����
                        string qrySeqs = informationsVec[4];  // �����qry����

                        outStruct.startMap[chromosome].push_back(refStart);  // ���������ʼλ��

                        vector<string> informationsVecTmp{&informationsVec[0], &informationsVec[0]+8};  // ��FORMAT�ĳ�GT����ߵ���Ϣȫ��ȥ��
                        outStruct.BaseInfoMap[chromosome][refStart] = join(informationsVecTmp, "\t") + "\tGT:SO:DP:ADP:NDP";  // �������informationsVecTmp��Ϣ GT:software:depth:averageDepth:normalDepth

                        // ��ȡ������Ϣ
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            {0, 0}, 
                            "all"
                        );

                        // ��������ref���Ⱥ�qry���ȡ�qryΪvector����Ϊ�ж����λ����
                        outStruct.refLenMap[chromosome].push_back(refLen);
                        outStruct.qryLenVecMap[chromosome].push_back(qryLenVec);

                        //���vector���ͷ��ڴ�
                        vector<int>().swap(qryLenVec);
                    }
                    // ����ַ���
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        return outStruct;
    }


    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index������base����
	 * @param vcfFileName        ��������vcf�ļ�
     * @param software           �����
     * @param sample_name        ��������������ȡ��������Ϣ
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
        softwareVcfStruct softvcfStructure;  // ���struct

        softvcfStructure.software = software;  // �����

        int sampleIdx = 0;  // ��¼sample��vcf�ļ��ĵڼ���

        cerr << "[" << __func__ << "::" << getTime() << "] " << "Building index: " << software << endl;  // print log

        vector<int> all_length_list;  // ѭ����vcf��Ϣ�ӵ�outStructure��

        gzFile gzfp = gzopen(vcfFileName.c_str(), "rb");  // �����ļ���

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << vcfFileName << "': No such file or directory." << endl;
            exit(1);
        }
        else
        {
            string informations;
            char line[1024]; // һ��ֻ��1024�ֽڵ�����

            while(gzgets(gzfp, line, 1024))
            {
                informations += line;

                if (informations.find("\n") != string::npos) // һ�н���
                {
                    if(informations.empty())
                    {
                        continue;
                    }

                    informations = strip(informations, '\n'); // ȥ�����з�

                    string::size_type findidx1 = informations.find("#");  // ע��������
                    if (findidx1 != string::npos)  // ע����
                    {
                        string::size_type findIdx2 = informations.find("#CHROM");  // #CHROMע��������
                        if (findIdx2 != string::npos)  // ����ǵĻ���sample_name����
                        {
                            vector<string> informationsVec = split(informations, "\t");  // vcfInfo���

                            // ��ȡsample��Ӧ������
                            if (software == "Paragraph")  // ��Ϊparagraph�ᱣ��ԭʼ���У����Դ�ԭʼ��֮����sample  mergeVcfStruct.colNum
                            {
                                sampleIdx = distance(informationsVec.begin(), find(informationsVec.begin()+mergeVcfStruct.colNum, informationsVec.end(), sample_name));  // ��ȡsample������λ��
                            }
                            else  // ����������ӵ�10�к��� FORMAT�ڵ�9�� 
                            {
                                sampleIdx = distance(informationsVec.begin(), find(informationsVec.begin()+9, informationsVec.end(), sample_name));  // ��ȡsample������λ��
                            }
                            
                            // ���û�ҵ��������˳�����
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

                            mergeVcfStruct.softwateSampleIdx[software] = sampleIdx; // �����Ӧ��sample�������
                        }
                    }
                    else  // ��ע����
                    {
                        // �ж�sample_name�ҵ�û��
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
                        string svType = informationsVec[2];  // �ж�BayesTyper�Ľ��Ϊduplication
                        string refSeq = informationsVec[3];
                        string qrySeqs = informationsVec[4];
                        string filter = informationsVec[6];
                        
                        // ��ȡλ��ĸ��Ƕ���Ϣ
                        vector<float> depthVec = get_depth(
                            informationsVec, 
                            sampleIdx
                        );
                        float depSum = accumulate(depthVec.begin(), depthVec.end(), 0);  // �������е�λ�������Ⱥ�
                        softvcfStructure.depthVec.push_back(depSum);  // �����Ӧ��vcf����б����
                        
                        vector<int> gtVec = get_gt(
                            informationsVec, 
                            sampleIdx
                        );
                        string gt = join(gtVec, "/");

                        // ���ݻ����ͽ��й���
                        if (gt == "0/0" || gt == "0" || gt == "." || gt == "./." || gtVec.size() == 0) //�����(0/0, .)��ʽ����ֱ�����������߷����˿��б�����������
                        {
                            // cerr << "[" << __func__ << "::" << getTime() << "] " 
                            //      << "Warning: skip -> " << informations << endl;
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // ����FILTER�ֶν��й���
                        if (filter != "PASS") //�������(PASS)��ֱ������������������
                        {
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                        
                        // ��������Ƿ�Խ��
                        int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
                        if (maxGtNum > split(qrySeqs, ",").size())  // �ȼ���Ƿ������Ƿ�Խ��
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " 
                                << "Error: number of genotyping and query sequences do not match -> " 
                                << informations << endl;
                            exit(1);
                        }

                        // ��ȡ�����͵ĳ�����Ϣ
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // ����������ʼ�ͷ�����Ϣ
                        softvcfStructure.startMap[chromosome].push_back(refStart);
                        softvcfStructure.gtVecMap[chromosome].push_back(gtVec);
                        softvcfStructure.depthVecMap[chromosome].push_back(depthVec);  // �����Ӧ��vcf��ȣ�ÿ����λ����ģ�

                        // ��������refLen��qryLenVec
                        softvcfStructure.refLenMap[chromosome].push_back(refLen);
                        softvcfStructure.qryLenVecMap[chromosome].push_back(qryLenVec);
                    }
                    // ����ַ���
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        // ���������ƽ�����Ƕȣ�����ͱ�׼��
        float aveDepth;
        float variance;
        float sd;
        tie(aveDepth, variance, sd) = cal_var_sd(softvcfStructure.depthVec);
        get<0>(mergeVcfStruct.depthMap[software]) = aveDepth;
        get<1>(mergeVcfStruct.depthMap[software]) = variance;
        get<2>(mergeVcfStruct.depthMap[software]) = sd;

        // ������Ľ���ӵ�ͼ��
        vcf_merge(mergeVcfStruct, softvcfStructure);

        return 0;
    }


    /**
	 * vcf_merge��recall����ӵ���ͼ�к���
	 *
     * @date 2023/3/7
	 * @param mergeVcfStruct     build_basefile_index������base����
	 * @param chromosome         Ⱦɫ���
     * @param trueRefStart       ��ʵ��ref��ʼλ��
     * @param software           ���
     * @param qryLenVec          ������͵ĵ����Ͷ�Ӧ����
     * @param gtVec              ����ķ���gt�б�
     * @param gt                 ��ǰѭ����gt
     * @param depth              λ���Ӧ�����
     * @param j                  ������������ѭ��������
     * @param k                  ���ֲ��ҷ���������ѭ��������
     * @param l                  ��ʵλ�㵥���͵�ѭ������
     * @param indexLeft          ���ֲ��ҷ���������
     * @param indexRight         ���ֲ��ҷ���������
     * 
     * 
     * @return tuple<int, int>   tuple<state, j>  0 -> ��ȷ��ӣ�-1 -> ��Ҫ������һ��ѭ��
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
        // �ȳ�ʼ��recall��ϣ��
        if (mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome].find(trueRefStart) == mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome].end())  // �����ʼλ�õ�һ�μ������ʼ��
        {
            mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart];
        }
        if (mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart].find(software) == mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart].end())  // ���Ⱦɫ���һ�μ������ʼ��
        {
            vector<float> depthVecTmp(qryLenVec.size());
            vector<int> gtVecTmp(qryLenVec.size());
            mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software] = make_tuple(depthVecTmp, gtVecTmp);
        }

        // gt��l����һ��ʱ0������һ�����ǣ�����һ��ѭ������ֹSNP�ж�ʱ�����
        if ((gt != 0 && l == 0) || (gt == 0 && l != 0))
        {
            return make_tuple(-1, j);
        }
        else
        {
            // ÿ�������͵������
            get<0>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;  // ��Ӧ������λ�����
            get<1>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;  // ��Ӧ������λ�û�����

            // ���Ķ��ֲ��ҵ����꣬��һ����λ�����ֱ����ӵ���refStart��
            indexLeft = k; 
            indexRight = k;

            if (j == 0 && gt == gtVec[1])  // ���������������һ���ģ�ֱ�Ӹ�ֵ��Ȼ��������һ��ѭ��
            {
                j++;
                get<0>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = depth;  // ��Ӧ������λ�����
                get<1>(mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome][trueRefStart][software])[j] = l;  // ��Ӧ������λ�û�����
            }
        }

        return make_tuple(0, j);
    }


    /**
	 * ������Ľ�����й��˲������������ļ�Ϊ���������ļ���
	 *
     * @date 2023/3/6
	 * @param mergeVcfStruct     build_basefile_index������base����
	 * @param softvcfStructure   build_softwarefile_index���������vcf����
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

        for (auto iter = softvcfStructure.startMap.begin(); iter != softvcfStructure.startMap.end(); iter++)  // ѭ�������vcf��Ⱦɫ�壩
        {
            string chromosome = iter->first;
            vector<int> startVec = softvcfStructure.startMap[chromosome];
            vector<int> refLenVec = softvcfStructure.refLenMap[chromosome];
            vector<vector<int> > qryLenVecVec = softvcfStructure.qryLenVecMap[chromosome];
            vector<vector<int> > gtVecVec = softvcfStructure.gtVecMap[chromosome];
            vector<vector<float> > depthVecVec = softvcfStructure.depthVecMap[chromosome];

            // �ȳ�ʼ��recall��ϣ��
            if (mergeVcfStruct.recallSoftwareGtDepVecMap.find(chromosome) == mergeVcfStruct.recallSoftwareGtDepVecMap.end())  // ���Ⱦɫ���һ�μ������ʼ��
            {
                mergeVcfStruct.recallSoftwareGtDepVecMap[chromosome];
            }

            if (find(chromosomeVec.begin(), chromosomeVec.end(), chromosome) == chromosomeVec.end())  // �ȿ�mergeVcfStructure����û�и�Ⱦɫ�壬û�еĻ������˳�����
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: chromosome is not present in the base vcf file -> " << chromosome << endl;
                exit(1);
            }
            else  // ������ڵĻ������ö��ֲ��ҷ�������ı��죬�ж������ǲ���ͬһ�����ǵĻ�push�����ǵĻ����һ���µ�
            {
                vector<int> mergeRefStartVec = mergeVcfStruct.startMap[chromosome];
                vector<int> mergeRefLenVec = mergeVcfStruct.refLenMap[chromosome];
                vector<vector<int> > mergeQryLenVecVec = mergeVcfStruct.qryLenVecMap[chromosome];  // vector<vector<qryLen>>

                for (int i = 0; i < startVec.size(); i++)  // ѭ�������vcf��λ�ã�
                {
                    int refStart = startVec[i];
                    int refLen = refLenVec[i];
                    vector<int> qryLenVec = qryLenVecVec[i];  // λ����Ͷ�Ӧ�ĵ�λ�������г���
                    vector<int> gtVec = gtVecVec[i];  // λ����Ͷ�Ӧ�ķ���
                    vector<float> depthVec = depthVecVec[i];  // λ����Ͷ�Ӧ�ĵ�λ�������

                    // ��¼���ֲ��ҷ��������������λ������ͬһ��refStart
                    int indexLeft = -1;
                    int indexRight = -1;

                    // �Զ����λ�������
                    for (size_t j = 0; j < qryLenVec.size(); j++)
                    {
                        int qryLen = qryLenVec[j];
                        int gt = gtVec[j];
                        float depth = depthVec[j];
                    
                        // ���ֲ��ҷ���200/10bp�ڵı���������
                        if (indexLeft == -1 && indexRight == -1)  // �µĵ�λ�����ٲ���
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
                            // true�������Ϣ
                            vector<int> mergeQryLenVec = mergeQryLenVecVec[k];  // �����λ����ĳ���
                            int trueRefLen = mergeRefLenVec[k];
                            int trueRefStart = mergeRefStartVec[k];

                            // ����ref��qry��ͬ�ĳ�������
                            vector<int> trueSeqLenVec;
                            trueSeqLenVec.push_back(trueRefLen);  // ����Ӳο�������ĳ���
                            trueSeqLenVec.insert(trueSeqLenVec.end(), mergeQryLenVec.begin(), mergeQryLenVec.end());  // �����

                            // �������Ϊ0���򱨴����˳����롣
                            if (mergeQryLenVec.size() == 0)
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " 
                                    << "Error: mergeQryLenVec.size() == 0 -> chromosome:" << chromosome 
                                    << " refStart:" << mergeRefStartVec[k] << endl;
                                exit(1);
                            }

                            for (size_t l = 0; l < trueSeqLenVec.size(); l++)  // һ��λ���ж����λ�����ʱ��trueSeqLenVec�洢�˸�����λ����ĳ��ȣ���˱���������λ��ı����Ƿ��merge��ߵ�һ�£�����0
                            {
                                int trueQryLen = trueSeqLenVec[l];

                                // ȱʧ
                                if (refLen >= 50 && qryLen < 50)
                                {
                                    if ((abs(refStart-trueRefStart)<=200) &&
                                        (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                        ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
                                    {
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
                                        goto stop;
                                    }
                                }
                                // ����
                                else if (refLen < 50 && qryLen >= 50)
                                {
                                    if ((abs(refStart-trueRefStart)<=200) &&
                                        (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                        ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
                                    {
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }

                                        goto stop;
                                    }
                                }
                                // �滻
                                else if (refLen >= 50 && qryLen >= 50)
                                {
                                    if ((abs(refStart-trueRefStart)<=200) &&
                                        (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                        ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                        ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
                                    {
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
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
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
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
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
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
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
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
                                        // ��ӵ���ͼ��
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

                                        // û���������һ��ѭ��
                                        if (recall_state == -1)
                                        {
                                            continue;
                                        }
                                        
                                        goto stop;
                                    }
                                }
                            }
                        }
                        stop:; // ����ҵ��ˣ����˳�Ƕ��ѭ���������õ����͵�ѭ����������һ��������
                    }
                }
            }
        }
        return 0;
    }


    /**
	 * ��ȡ����ĳ�����Ϣ
	 *
     * @param vcfInfo           vcfInfo
     * @param gtVec             �������б�
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

        // ��ȡref�͵����͵ĳ���
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
            if (gtVec[i] == 0)  // �����������0������
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
	 * ��ȡ�����Ͷ�Ӧ�ĳ�����Ϣ.
	 *
     * @param svType                         ��������
	 * @param refSeq                         ref����Ϣ
     * @param qrySeqs                        qry����Ϣ
     * @param gtVec                          λ��ķ�����Ϣ
     * @param lenType                        ֻȡ�����Ͷ�Ӧ�ĳ��Ȼ������еĳ���(hap/all)
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
        // ���ģʽ�Ƿ���ȷ������ȷ�˳�����
        if (lenType != "hap" && lenType != "all")
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "Error: lenType -> " 
                << lenType << endl;
            exit(1);
        }
        
        int refLen = refSeq.size();  // ref���г���
        vector<string> qrySeqVec = split(qrySeqs, ",");  // qry�������б�
        
        // ��������Ƿ�Խ�磬���ֻҪhap�ĳ���ʱ���ټ��
        if (lenType == "hap")
        {
            int maxGtNum = *max_element(gtVec.begin(), gtVec.end());
            if (maxGtNum > qrySeqVec.size())  // �ȼ���Ƿ������Ƿ�Խ��
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: number of genotyping and query sequences do not match -> " 
                    << qrySeqs << endl;
                exit(1);
            }
        }

        // ����ref��qry��ͬ�ĳ�������
        vector<int> seqLenVec;
        seqLenVec.push_back(refLen);

        // ������λ�����б�
        for (size_t i = 0; i < qrySeqVec.size(); i++)
        {
            string qrySeq = qrySeqVec[i];

            // ��ʱ�洢����
            int refLenTmp = refLen;
            int qryLenTmp = qrySeq.size();

            string::size_type idxSvsize = qrySeq.find(">");
            string::size_type idxIns = qrySeq.find("<INS");  // <INS:SVSIZE=90:BREAKPOINT1>
            string::size_type idxDel = qrySeq.find("<DEL");  // <DEL:SVSIZE=1720:BREAKPOINT>
            string::size_type idxDup = qrySeq.find("<DUP");  // <DUP:SVSIZE=10001:COVERAGE>

            // �ж�BayesTyper�Ľ��Ϊduplication
            string::size_type idxSvTypeDup = svType.find("Duplication");

            if (idxSvsize != string::npos)  // GraphTyper2�Ľ��
            {
                if (idxIns != string::npos && idxDel == string::npos) // �ַ����а�������
                {
                    refLenTmp = refLen;
                    qryLenTmp = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                }
                else if (idxIns == string::npos && idxDel != string::npos) // �ַ����а���ȱʧ
                {
                    refLenTmp = stoi(split(split(qrySeq, "=")[1], ":")[0]);
                    qryLenTmp = refLen;
                }
                else if (idxDup != string::npos) // �ַ����а����ظ����ֶΣ�GraphTyper2��<DUP:SVSIZE=2806:BREAKPOINT1>
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

                seqLenVec[0] = refLenTmp; // ����ref�ĳ���
            }

            // �ж�BayesTyper�Ľ��Ϊduplication����ref_lenΪ1-2��qry_seqȴ�ܳ�ʱ
            if (idxSvTypeDup != string::npos && refLenTmp <= 2) // BayesTyper���duplication��ɲ��룬�������¼��㳤��
            {
                refLenTmp = qryLenTmp;
                qryLenTmp *= 2;
            }

            // ���qry�ĳ���
            if (seqLenVec.size() == (i + 1))  // �жϵ�����ȷʵ���Լ���λ���ϣ�
            {
                seqLenVec.push_back(qryLenTmp);  // ���qry�ĳ���
            }
            else  // �������б��Ȳ�����ʱ����
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "Error: wrong index for haplotype length -> " 
                    << qrySeqs 
                    << endl;
                exit(1);
            }
        }

        // �һ����Ͷ�Ӧ��qry����
        vector<int> qryLenVec;
        if (lenType == "hap")  // ֻȡ�����͵����г���
        {
            if (gtVec.size() == 1)  // gtֻ��һ��������Ϊ���ϵı��죬�������������һ���ĳ���
            {
                if (gtVec[0] == 0)  // ���������Ϊ0��������λ��
                {
                    // ���� '0/0'
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
        else  // ȡ���еĳ���
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
	 * @param gtAllVec        λ�����е�gt��vector
	 * @param softwareVec     λ�����е�software��vector
	 * @param depthNorVec     λ�����б�׼�����depth��vector
     * @param softwareNum     �û�������������
     * @param vcfInfo         λ��ı�����Ϣ
     * @param mode            EVG���е�ģʽ
     * @param mergeVcfStruct  ��ͼ
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
        // �������
        vector<string> bestSoftwareVec{};  // ����Ƶ����ߵ����vector
        string selectSoftware;  // ����sv�ĳ��Ƚ���ɸѡ������vector
        int svLength{};
        int gtFrequency = 0;

        // ��׼����������ӽ�0������
        float selectDepth = 10000.0;
        int selectDepthIdx = -1;
        for (int i = 0; i < depthNorVec.size(); i++)
        {
            if (abs(depthNorVec[i]) < selectDepth)
            {
                selectDepth = abs(depthNorVec[i]);  // �������
                selectDepthIdx = i;  // ��������
            }
        }

        // ����ÿһ�ֻ����ͳ��ֵ�Ƶ��
        // ����Ƶ������gt
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

        // ��Ƶ������gt��Ӧ����ҳ���
        for (size_t i = 0; i < gtAllVec.size(); i++)
        {
            if (count(gtAllVec.begin(), gtAllVec.end(), gtAllVec[i]) == gtFrequency)
            {
                bestSoftwareVec.push_back(softwareVec[i]);
            }
        }

        // ���˹���
        if (-49 <= svLength && svLength <= 49)  // snp+indel
        {
            if (mode == "specific")  // Ϊ��ϵ�����vcf�ļ�
            {
                if (gtFrequency < 2) // ֧�ֵ��������С��2��ѡȡ�����ӽ���0�Ļ�����
                {
                    selectSoftware = softwareVec[selectDepthIdx];  // ѡȡ���ǶȽӽ���0�����
                }
                else // ����������������֧�ֵĻ�����ѡ��û����͡�
                {
                    selectSoftware = bestSoftwareVec[0]; // ѡ�����Ƶ����ߵķ��ͽ��
                }
            }
            else  // Ϊ��ͼ�õ�vcf
            {
                if (gtFrequency < 2) // ����������������ϵ����֧�֣�����������λ��
                {
                    selectSoftware.clear();
                }
                else // ����������������֧�ֵĻ�����ѡ��û����͡�
                {
                    selectSoftware = bestSoftwareVec[0]; // ѡ�����Ƶ����ߵķ��ͽ��
                }
            }
        }
        else  // Insertion and Deletion
        {
            vector<string> softwareSortVec = {"BayesTyper", "Paragraph", "VG-MAP", "VG-Giraffe", "GraphTyper2", "GraphAligner", "PanGenie"};
            for (auto it : softwareSortVec)
            {
                // ����50bp�ı����������it֧�֣���ֱ��ѡ��it�Ľ��
                if (find(softwareVec.begin(), softwareVec.end(), it) != softwareVec.end())
                {
                    selectSoftware = it;
                    break;
                }
                else // ����ѡȡ���ǶȽӽ�0�����
                {
                    selectSoftware = softwareVec[selectDepthIdx];  // ѡȡ���ǶȽӽ���0�����
                }
            }
        }

        return selectSoftware;
    }



    /*
        vcf����
        mergeVcfStruct -> ����ϲ����struct
    */

   /**
	 * Filter and merge result according to the depth and number of variants.
	 *
	 * @param mergeVcfStruct     ����ϲ����struct
	 * @param softwareNum        �û��ύ���������
	 * @param mode               EVG���е�ģʽ
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

        // ������
        map<string, map<int, string > > outMap;  // outMap<chromosome, map<refStart, vcfInfo> >

        // ����recall��ϣ��
        for (auto it1 : mergeVcfStruct.recallSoftwareGtDepVecMap)  // map<chr, map<start, map<software, tuple<vector<depth>, vector<gt> > > > >
        {
            string chromosome = it1.first;
            
            for (auto it2 : it1.second) // map<start, map<software, tuple<vector<depth>, vector<gt> > > >
            {
                int refStart = it2.first;

                // ��λ���vcfInfo
                string informations = mergeVcfStruct.BaseInfoMap[chromosome][refStart];

                // ����filter_rule�Ĳ���
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

                    // ��ȡ���ƽ���ĸ��Ƕȣ��ȼ��map����û�У����û���򱨴�
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
                    float depthSumTmp = accumulate(depthVecTmp.begin(), depthVecTmp.end(), 0);  // �������е�λ�������Ⱥ�
                    depthSumVec.push_back(depthSumTmp);
                    
                    // �Ը��Ƕȱ�׼����z-score��׼����
                    float depthNor = (depthSumTmp - meanDepth) / (float)sd;
                    
                    gtAllVec.push_back(gtVecTmp);
                    softwareVec.push_back(softwareTmp);
                    depthNorVec.push_back(depthNor);
                }

                // ���ˣ�ѡ�����п��ܵĵ�����
                string selectSoftware = filter_rule(
                    gtAllVec, 
                    softwareVec, 
                    depthNorVec, 
                    informations, 
                    softwareNum, 
                    mode, 
                    mergeVcfStruct
                );

                if (selectSoftware.length() == 0) // ������ص��ǿյ��ַ����������λ�㲻���ţ�������λ��
                {
                    continue;
                }

                // ѡ�����������������ͺ������Ϣ
                int softwareIdx = distance(softwareVec.begin(), find(softwareVec.begin(), softwareVec.end(), selectSoftware));  // ��ȡsoftware������λ��
                vector<int> selectGtVec = gtAllVec[softwareIdx];
                float selectDepthSumVec = depthSumVec[softwareIdx];
                float selectDepthNor = depthNorVec[softwareIdx];
                vector<float> selectDepthVec = depthVecVec[softwareIdx];

                // // ������
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
	 * ���㷽��ͱ�׼��
	 *
	 * @param data     ����λ����ȵ��б�
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
        ������
        mergeVcfStruct -> ����ϲ����struct
        outMap -> vcf_merge_filter������   outMap[chromosome][refStart][refLen][qryLen]
        prefix -> ����ļ���ǰ׺
    */

    /**
	 * ������
	 *
     * @param headInfo           vcfע����
	 * @param outMap             vcf_merge_filter������   outMap[chromosome][refStart]
	 * @param outputFileName     ����ļ���
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

        // ������Ҫ�����string
        string outTxt = headInfo;

        for (auto it1 : outMap)  // outMap<chromosome, map<refStart, vcfInfo> >
        {
            for (auto it2: it1.second)  // map<refStart, vcfInfo>
            {
                string vcfInfo = it2.second;
                outTxt += vcfInfo + "\n";
            }
        }

        if (outputFileName.size() == 0)  // û��ָ������ļ�����ӡ����Ļ
        {
            cout << outTxt << endl;
        }
        else if (outputFileName.find(".gz") != string::npos)  // ���Ϊѹ���ļ�
        {
            // ����ļ���
            gzFile gzfp = gzopen(outputFileName.c_str(), "wb");

            // ���ļ�
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
            // �ͷ��ڴ棬�ر��ļ�
            gzclose(gzfp);
        }
        else
        {
            // ����ļ���
            ofstream outputFile;
            outputFile.open(outputFileName, ios::out);

            // ���ļ�
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

            // �ر��ļ�
            outputFile.close();
        }

        return 0;
    }
}  // namespace MERGE

#endif