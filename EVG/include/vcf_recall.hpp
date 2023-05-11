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
    // ȫ�ֱ���
    // rocKey  
    string rocKey;  // (FORMAT.<name>/INFO.<name>) [FORMAT.DP]
    // rocType�ж�����recall������roc������ֻ��genotype����roc  recall/genotype [genotype]
    string rocType;
    

    struct vcfStructure
    {
        map<string,vector<int> > refStartVecMap;  // map<chr,vector<refStart> >
        
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > chrStartLenInfoGtVecMap;  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen, gtVec> > >
        
        vector<int> allLengthList;  // ���б���ĳ���
    };


    /**
	 * roc�����޸�.
	 *
	 * @param rocKeyTmp        main�����rocKey
     * @param rocTypeTmp       main�����rocTypeTmp
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


    // �Ա��쳤�Ƚ���ͳ�Ƶ�Bool����
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
	 * ��vcf�ļ��������򲢹�������
	 *
	 * @param genotype_filename    �����ļ�
     * 
     * 
     * @return outStructure        vcfStructure
	**/
    vcfStructure build_index(
        string genotype_filename
    )
    {
        vcfStructure outStructure;

        // ���vcf�ļ��Ƿ�����
        check_vcf_sort(genotype_filename);

        // �����ļ���
        gzFile gzfp = gzopen(genotype_filename.c_str(), "rb");

        // ѭ����genotype��ӵ�genotype vector��
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

                    // ����ע����
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo���
                        string chromosome = informationsVec[0];  // Ⱦɫ����Ϣ
                        int refStart = stoi(informationsVec[1]);  // �������ʼλ��
                        string svType = informationsVec[2];  // �ж�BayesTyper�Ľ��Ϊduplication
                        string refSeq = informationsVec[3];  // �����ref����
                        string qrySeqs = informationsVec[4];  // �����qry����

                        vector<int> gtVec = get_gt(
                            informationsVec
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

                        // ��ȡref�͵����͵ĳ���
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // ��ȡ����ĳ���
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // �������б���ĳ���
                        outStructure.allLengthList.push_back(svLength);

                        // ������Ϣ
                        outStructure.refStartVecMap[chromosome].push_back(refStart);

                        // �ȳ�ʼ����ϣ��
                        if (outStructure.chrStartLenInfoGtVecMap.find(chromosome) == outStructure.chrStartLenInfoGtVecMap.end())
                        {
                            outStructure.chrStartLenInfoGtVecMap[chromosome];
                        }
                        outStructure.chrStartLenInfoGtVecMap[chromosome][refStart] = make_tuple(refLen, qryLenVec, informations, svLength, gtVec);
                    }

                    // ����ַ���
                    informations.clear();
                    string().swap(informations);
                }
            }
        }
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        return outStructure;
    }


    /**
	 * ��gramtools�Ľ������ת��
	 *
	 * @param evaluateFilename    �����ļ�
     * 
     * 
     * @return outFileName        ����ļ�
	**/
    string gramtools_convert(
        string evaluateFilename
    )
    {
        string informations;
        string outInformationsTitle;
        string outInformations;

        // ����ļ��Ƿ�����
        check_vcf_sort(
            evaluateFilename
        );

        // ������
        // ����ļ���
        vector<string> prefixVec = split(evaluateFilename, "/");
        string outFileName = "convert." + prefixVec[prefixVec.size()-1];
        if (outFileName.find(".gz") == string::npos && outFileName.find(".GZ") == string::npos)
        {
            outFileName += ".gz";
        }
        
        gzFile gzfp1 = gzopen(outFileName.c_str(), "wb");

        // ���ļ�
        if(!gzfp1)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << outFileName 
                << "': No such file or directory." 
                << endl;
            exit(1);
        }


        // �����ļ���
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
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

                    // ����ע����
                    if (informations.find("#") == string::npos)
                    {
                        // ����ע����
                        if (outInformationsTitle.size() > 0)
                        {
                            gzwrite(gzfp1, outInformationsTitle.c_str(), outInformationsTitle.length());

                            // ���ע�����ַ���
                            outInformationsTitle.clear();
                            string().swap(outInformationsTitle);
                        }
                        // ÿ10Mbд��һ���ļ�
                        if (outInformations.size() > 10000000)
                        {
                            gzwrite(gzfp1, outInformations.c_str(), outInformations.length());

                            // ����ַ���
                            outInformations.clear();
                            string().swap(outInformations);
                        }

                        // ���ַ������зָ���浽vector��
                        vector<string> informationsVec = split(informations, "\t");

                        // ��ȡinformationsVector�л����ͺ�λ�㸲�Ƕ���Ϣ
                        vector<string> gtVector = split(strip(informationsVec[informationsVec.size()-1], '\n'), ":");
                        vector<string> formatVector = split(strip(informationsVec[8], '\n'), ":");
                        string gt = gtVector[0];

                        // ��FORMAT�ֶ���COV��λ��
                        vector<string>::iterator covItera = find(formatVector.begin(), formatVector.end(), "COV");
                        int covIndex = 0;
                        if (covItera != formatVector.end()) // FORMAT����COV
                        {
                            covIndex = distance(formatVector.begin(), covItera);
                        }
                        else // û�еĻ�ֱ��ת����������һ��ѭ��
                        {
                            outInformations += informations + "\n";
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // FILTER�ֶ�ǰ�ߵ����ݲ��䣬ֱ�Ӽӵ�outInformations��
                        for (int i = 0; i < informationsVec.size()-1; i++)
                        {
                            if (i == 6) // FILTER�и�ΪPASS��Դ�ļ��е�Ϊ.
                            {
                                outInformations += "PASS\t";
                            }
                            else
                            {
                                outInformations += informationsVec[i] + "\t";
                            }
                        }

                        string siteCoverage = gtVector[covIndex];
                        

                        // ���gt=0�Ļ������λ��û�б��죬��gt��Ϊ0/0
                        if (gt == "0")
                        {
                            outInformations += "\t0/0";
                        }
                        // gt��Ϊ0��ʱ��
                        else
                        {
                            // ���λ�㸲�Ƕ�����0,�ֶΣ������λ��Ϊ���͵�ͻ��λ��
                            if (siteCoverage.find("0,") < siteCoverage.length())
                            {
                                outInformations += "\t1/1";
                            }
                            // ���λ�㸲�Ƕ�û��0,�ֶΣ������λ��Ϊ�Ӻϵ�ͻ��λ��
                            else
                            {
                                outInformations += "\t0/1";
                            }
                        }

                        // ��informations�����ֶμ���
                        for (int i = 1; i < gtVector.size(); i++)
                        {
                            outInformations += ":" + gtVector[i];
                        }

                        outInformations += "\n";
                    }
                    // ����ע������Ϣ
                    else
                    {
                        outInformationsTitle += informations + "\n";
                    }
                    
                    // ����ַ���
                    informations.clear();
                    string().swap(informations);
                }
            }
        }

        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        // ���д��һ���ļ�
        gzwrite(gzfp1, outInformations.c_str(), outInformations.length());
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp1);
        
        return outFileName;
    }


    /**
	 * genotype��������
	 *
	 * @param evaluateFilename            �����ļ�����������
     * @param chrStartLenInfoGtVecMap     �漯���е�vcf��Ϣ
     * @param all_length_list             �漯����vcf����
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
        
        // ���쳤�ȷ�Χ
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
        map<float, vector<int> > rocCallMap; // map<DP/GQ, vector<length>> ����ҵ����е�
        map<float, vector<int> > rocRecallMap; // map<DP/GQ, vector<length>> ����ҵ���ȷ��
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

        // ����ļ��Ƿ�����
        check_vcf_sort(evaluateFilename);
        
        // ���������ȷ�Ľ��
        string true_txt;
        gzFile trueFile = gzopen("genotype.true.vcf.gz", "wb");
        // ���ļ�
        if(!trueFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ������ʹ���Ľ��
        string genotypeMisTxt;
        gzFile genotypeMisFile = gzopen("genotype.err.vcf.gz", "wb");
        // ���ļ�
        if(!genotypeMisFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.err.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ����miscall�Ľ��
        string misCallTxt;
        gzFile misCallFile = gzopen("miscall.vcf.gz", "wb");
        // ���ļ�
        if(!misCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'miscall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // û���ҵ��ı���
        string failCallTxt;
        gzFile failCallFile = gzopen("failcall.vcf.gz", "wb");
        // ���ļ�
        if(!failCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'failcall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ����vector
        vector<int> genotype_length_list;
        vector<int> misgenotype_length_list;
        vector<int> miscall_length_list;

        // �����ļ���
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
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

                    // ����ע����
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");  // vcfInfo���
                        string chromosome = informationsVec[0];
                        int refStart = stoi(informationsVec[1]);  // �������ʼλ��
                        string svType = informationsVec[2];  // �ж�BayesTyper�Ľ��Ϊduplication
                        string refSeq = informationsVec[3];  // �����ref����
                        string qrySeqs = informationsVec[4];  // �����qry����
                        string filter = informationsVec[6];

                        // ����FILTER�ֶν��й���
                        if (filter != "PASS") // �������(PASS)��ֱ������
                        {
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // ��ȡ��������Ϣ
                        vector<int> gtVec = get_gt(
                            informationsVec
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

                        /* --------------------------------- roc ------------------------------------ */
                        vector<string> lastColVec = split(informationsVec[informationsVec.size()-1], ":");
                        float rocNum;
                        if (rocColNum == 8) // FORMAT�ֶ�
                        {
                            vector<string> rocInfVec = split(strip(informationsVec[rocColNum], '\n'), ":");
                            // rocNum���±�
                            vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec[1]);
                            if (rocItera == rocInfVec.end()) // �����ֶ�����û�ж�Ӧ��rocNum��Ϣ
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
                        }
                        else // INFO�ֶ�
                        {
                            smatch patternResult; // ������ʽ�Ľ��
                            string rocInfString = strip(informationsVec[rocColNum], '\n');
                            regex pattern(rocKeyVec[1]  + "=(\\d+)");
                            string::const_iterator iterStart = rocInfString.begin();
                            string::const_iterator iterEnd = rocInfString.end();
                            regex_search(iterStart, iterEnd, patternResult, pattern);
                            string rocNumString = patternResult[1];

                            // �����û�н����û�еĻ�����
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

                        // ��ȡref�͵����͵ĳ���
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // ��ȡ����ĳ���
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // ��ӳ��ȵ�rocAllMap��
                        rocCallMap[rocNum].push_back(svLength);

                        // ���л����������
                        // ���Ⱦɫ���Ƿ���ڣ������ھ���miscall
                        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > >::iterator iter1 = chrStartLenInfoGtVecMap.find(chromosome);
                        if (iter1 == chrStartLenInfoGtVecMap.end())
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Warning: " << chromosome << " not in chrStartLenInfoGtVecMap.\n";
                            miscall_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // ��¼genotype�Ľ��
                        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype
                        int trueRefLen;
                        string trueVcfInfo;
                        int trueSvLen;

                        // �ڶ�ӦȾɫ���map�в�����ֵ��refLen qryLenVec
                        map<int, tuple<int, vector<int>, string, int, vector<int> > >::iterator iter2 = iter1->second.find(refStart);
                        if (iter2 != iter1->second.end())  // ����ҵ��ˣ��������жϷ����Ƿ���ȷ
                        {
                            // true������Ϣ
                            trueRefLen = get<0>(iter2->second);
                            vector<int> trueQryLenVec = get<1>(iter2->second);
                            trueVcfInfo = get<2>(iter2->second);
                            trueSvLen = get<3>(iter2->second);
                            vector<int> trueGtVec = get<4>(iter2->second);

                            // �������Ϊ0���򱨴����˳����롣
                            if (trueQryLenVec.size() == 0)
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] " 
                                    << "Error: trueQryLenVec.size() == 0 -> chromosome:" << chromosome 
                                    << " refStart:" << refStart << endl;
                                exit(1);
                            }

                            // ������λ����
                            int trueQryLenVecIdx = -1;  // ��¼��ֵ���ù��ĵ����ͣ���ֹ�ظ�����
                            for (size_t i = 0; i < qryLenVec.size(); i++)
                            {
                                int qryLen = qryLenVec[i];
                                int callGt = gtVec[i];

                                for (size_t j = 0; j < trueQryLenVec.size(); j++)
                                {
                                    // �ù��ĵ���������
                                    if (j == trueQryLenVecIdx)
                                    {
                                        continue;
                                    }

                                    int trueQryLen = trueQryLenVec[j];
                                    int trueGt = trueGtVec[j];

                                    // callGt��trueGt����һ����0������һ�����ǣ�����һ��ѭ������ֹSNP�ж�ʱ�����
                                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0))
                                    {
                                        continue;;
                                    }

                                    // ȱʧ
                                    if (refLen >= 50 && qryLen < 50)
                                    {
                                        if ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;

                                            goto stop;
                                        }
                                    }
                                    // ����
                                    else if (refLen < 50 && qryLen >= 50)
                                    {
                                        if ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25)
                                        {
                                            genotypeTrueNum++;
                                            trueQryLenVecIdx = j;
                                            
                                            goto stop;
                                        }
                                    }
                                    // �滻
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
                                stop:; // ����ҵ��ˣ����˳�Ƕ��ѭ���������õ����͵�ѭ����������һ��������
                            }
                        }

                        // �жϷ��͵Ľ��
                        if (genotypeTrueNum == 0)  // �ҵ��ı��첻���漯��
                        {
                            miscall_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                        else
                        {
                            // ������ȷ
                            if (genotypeTrueNum >= gtVec.size())
                            {
                                genotype_length_list.push_back(trueSvLen);
                                true_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                            "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                                // ��ӳ��ȵ�rocTrueMap��
                                rocRecallMap[rocNum].push_back(svLength);
                            }
                            // ���ʹ���
                            else
                            {
                                misgenotype_length_list.push_back(trueSvLen);
                                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                                "true_length\t" + to_string(trueSvLen) + "\t" + trueVcfInfo + "\n";

                                // ���ж�����recall������roc�������������ֻ��genotype����roc
                                if (rocType == "recall")
                                {
                                    // ��ӳ��ȵ�rocTrueMap��
                                    rocRecallMap[rocNum].push_back(svLength);
                                }
                                
                            }
                            // ɾ���Ѿ����ù��ı���
                            chrStartLenInfoGtVecMap[chromosome].erase(refStart);
                        }
                    }
                    // ����ַ���
                    informations.clear();
                    string().swap(informations);

                    // ������
                    if (true_txt.size() > 10000000 || 
                        genotypeMisTxt.size() > 10000000 || 
                        misCallTxt.size() > 10000000) // ÿ10Mbд��һ��
                    {
                        gzwrite(trueFile, true_txt.c_str(), true_txt.length()); // ������ȷ
                        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // ���ʹ���
                        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // �����漯�еı���

                        // ����ַ���
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
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        // ������
        gzwrite(trueFile, true_txt.c_str(), true_txt.length()); // ������ȷ

        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // ���ʹ���

        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // �����漯�еı���

        // δ�ҵ�����ʵ����
        for (auto it1 : chrStartLenInfoGtVecMap)  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        {
            for (auto it2 : it1.second)  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            {
                failCallTxt += get<3>(it2.second) + "\n";

                if (failCallTxt.size() > 10000000) // ÿ10Mbд��һ��
                {
                    gzwrite(failCallFile, failCallTxt.c_str(), failCallTxt.length()); // ������ȷ

                    // ����ַ���
                    failCallTxt.clear();
                    string().swap(failCallTxt);
                }
            }
        }
        gzwrite(failCallFile, failCallTxt.c_str(), failCallTxt.length());

        // �ر��ļ�
        gzclose(trueFile);
        gzclose(genotypeMisFile);
        gzclose(misCallFile);
        gzclose(failCallFile);

        int sv_genotype_recall = genotype_length_list.size();
        int sv_misgenotype_recall = misgenotype_length_list.size();
        int sv_mis_call = miscall_length_list.size();
        int sv_all = all_length_list.size();

        // ������
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

        // ���Ƚ�����
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

        // �ͷ��ڴ�
        outFile.close();

        roc_calculate(rocCallMap, rocRecallMap, all_length_list);

        return 0;
    }


    /**
	 * recall��������
	 *
	 * @param evaluateFilename            �����ļ�����������
     * @param trueVcfStructure            build_index()������
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

        // ���쳤�ȷ�Χ
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
        map<float, vector<int>> rocCallMap; // map<DP/GQ, vector<length>> ����ҵ����е�
        map<float, vector<int>> rocRecallMap; // map<DP/GQ, vector<length>> ����ҵ���ȷ��
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

        // ���vcf�ļ��Ƿ�����
        check_vcf_sort(evaluateFilename);

        // ������ȷ�Ľ��
        string true_txt{};
        gzFile trueFile = gzopen("recall.true.vcf.gz", "wb");
        // ���ļ�
        if(!trueFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'recall.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ���������ȷ�Ľ��
        string true_Gt_txt{};
        gzFile trueGtFile = gzopen("genotype.true.vcf.gz", "wb");
        // ���ļ�
        if(!trueGtFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.true.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ������ʹ���Ľ��
        string genotypeMisTxt{};
        gzFile genotypeMisFile = gzopen("genotype.err.vcf.gz", "wb");
        // ���ļ�
        if(!genotypeMisFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'genotype.err.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ����miscall�Ľ��
        string misCallTxt{};
        gzFile misCallFile = gzopen("miscall.vcf.gz", "wb");
        // ���ļ�
        if(!misCallFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'miscall.vcf.gz'" 
                << ": No such file or directory." 
                << endl;
            exit(1);
        }

        // ����vector
        vector<int> genotype_length_list;
        vector<int> call_length_list;
        vector<int> recall_length_list;

        // ����ע����
        string start_chromosome; // �����ж��ǲ����µ�Ⱦɫ�壬�ǵĻ������¹���vector��

        // ��¼�Ѿ����ù��ı��죬��ֹ�ظ�����
        map<string,vector<int>> selectChrStartMap;

        vector<int> start_vector;
        vector<int> ref_len_vector;
        vector<int> qry_len_vector;
        vector<string> gt_vector;
        vector<string> vcf_inf_vector;

        // ���ֲ��ҷ���������
        int leftIdxTmp = 0;
        int rightIdxTmp = 0;

        // �����ļ���
        gzFile gzfp = gzopen(evaluateFilename.c_str(), "rb");

        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << evaluateFilename << "': No such file or directory." << endl;
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

                    // ����ע����
                    if (informations.find("#") == string::npos)
                    {
                        vector<string> informationsVec = split(informations, "\t");

                        string chromosome = informationsVec[0];
                        int refStart = stoi(informationsVec[1]);  // �������ʼλ��
                        string svType = informationsVec[2];  // �ж�BayesTyper�Ľ��Ϊduplication
                        string refSeq = informationsVec[3];  // �����ref����
                        string qrySeqs = informationsVec[4];  // �����qry����
                        string filter = informationsVec[6];

                        // ����FILTER�ֶν��й���
                        if (filter != "PASS") // �������(PASS)��ֱ������
                        {
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }

                        // ��ȡ��������Ϣ
                        vector<int> gtVec = get_gt(
                            informationsVec
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


                        /* --------------------------------- roc ------------------------------------ */
                        vector<string> lastColVec = split(informationsVec[informationsVec.size()-1], ":");
                        float rocNum;
                        if (rocColNum == 8) // FORMAT�ֶ�
                        {
                            vector<string> rocInfVec = split(strip(informationsVec[rocColNum], '\n'), ":");
                            // rocNum���±�
                            vector<string>::iterator rocItera = find(rocInfVec.begin(), rocInfVec.end(), rocKeyVec[1]);
                            if (rocItera == rocInfVec.end()) // �����ֶ�����û�ж�Ӧ��rocNum��Ϣ
                            {
                                cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Error: " << rocKeyVec[1] << " not in " << rocKeyVec[0] << " column -> " 
                                    << " (" << strip(informations, '\n') << ")"
                                    << endl;
                                exit(1);
                            }
                            rocNum = stof(lastColVec[distance(rocInfVec.begin(), rocItera)]);
                        }
                        else // INFO�ֶ�
                        {
                            smatch patternResult; // ������ʽ�Ľ��
                            string rocInfString = strip(informationsVec[rocColNum], '\n');
                            regex pattern(rocKeyVec[1]  + "=(\\d+)");
                            string::const_iterator iterStart = rocInfString.begin();
                            string::const_iterator iterEnd = rocInfString.end();
                            regex_search(iterStart, iterEnd, patternResult, pattern);
                            string rocNumString = patternResult[1];

                            // �����û�н����û�еĻ�����
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


                        // ��ȡref�͵����͵ĳ���
                        int refLen;
                        vector<int> qryLenVec;
                        tie(refLen, qryLenVec) = get_hap_len(
                            svType, 
                            refSeq, 
                            qrySeqs, 
                            gtVec, 
                            "hap"
                        );

                        // ��ȡ����ĳ���
                        int svLength = sv_length_select(
                            refLen, 
                            qryLenVec, 
                            gtVec
                        );

                        // ��ӳ��ȵ�rocAllMap��
                        rocCallMap[rocNum].push_back(svLength);


                        // recall
                        // ���Ⱦɫ���Ƿ���ڣ������ھ���miscall
                        if (trueVcfStructure.chrStartLenInfoGtVecMap.find(chromosome) == trueVcfStructure.chrStartLenInfoGtVecMap.end())
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] "
                                    << "Warning: " << chromosome << " not in chrStartLenInfoGtVecMap.\n";
                            call_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                            // ����ַ���
                            informations.clear();
                            string().swap(informations);
                            continue;
                        }
                    

                        if (start_chromosome != chromosome)  // ���ֲ��ҵ������������㣬���ڼӿ��ѯ�ٶ�
                        {
                            // ������������
                            leftIdxTmp = 0;
                            rightIdxTmp = 0;
                        }

                        // ��¼���ֲ��ҷ��������������λ������ͬһ��refStart
                        int indexLeft = -1;
                        int indexRight = -1;

                        int genotypeTrueNum = 0;  // 0->misCall >0&&<gtVec.size()->recall  >=gtVec.size()->trueGenotype

                        // ��¼recall��  trueRefStart, trueRefLen, trueSvLen��truevcfInfo
                        int trueRefStart;
                        int trueRefLen;
                        int trueSvLen;
                        string truevcfInfo;

                        // ��¼�漯�б��ù��ĵ�����
                        int trueQryLenVecIdx = -1;

                        // �������λ����
                        for (size_t i = 0; i < qryLenVec.size(); i++)
                        {
                            int qryLen = qryLenVec[i];
                            int callGt = gtVec[i];

                            // ���ֲ��ҷ���200/10bp�ڵı���������
                            if (indexLeft == -1 && indexRight == -1)  // �µĵ�λ�����ٲ���
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

                            // �������ֲ��ҷ�������
                            for (int j = indexLeft; j <= indexRight; j++)
                            {
                                // true�������Ϣ
                                trueRefStart = trueVcfStructure.refStartVecMap[chromosome][j];

                                // �ȼ�������û�б��ù��������ɾ���˾���һ������
                                if (trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].find(trueRefStart) == trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].end())
                                {
                                    continue;
                                }
                                
                                trueRefLen = get<0>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                vector<int> trueQryLenVec = get<1>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                truevcfInfo = get<2>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                trueSvLen = get<3>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);
                                vector<int> trueGtVec = get<4>(trueVcfStructure.chrStartLenInfoGtVecMap[chromosome][trueRefStart]);

                                for (size_t k = 0; k < trueQryLenVec.size(); k++)  // һ��λ���ж����λ�����ʱ��trueSeqLenVec�洢�˸�����λ����ĳ��ȣ���˱���������λ��ı����Ƿ��merge��ߵ�һ�£�����0
                                {
                                    // �ù��ĵ���������
                                    if (k == trueQryLenVecIdx)
                                    {
                                        continue;
                                    }

                                    int trueQryLen = trueQryLenVec[k];
                                    int trueGt = trueGtVec[k];

                                    // callGt��trueGt����һ����0������һ�����ǣ�����һ��ѭ������ֹSNP�ж�ʱ�����
                                    if ((callGt != 0 && trueGt == 0) || (callGt == 0 && trueGt != 0))
                                    {
                                        continue;;
                                    }

                                    // ȱʧ
                                    if (refLen >= 50 && qryLen < 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
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
                                    // ����
                                    else if (refLen < 50 && qryLen >= 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // �滻
                                    else if (refLen >= 50 && qryLen >= 50)
                                    {
                                        if ((abs(refStart-trueRefStart)<=200) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=200) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // snp
                                    else if (refLen == 1 && qryLen == 1)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    // indel-del
                                    else if (refLen >= 3 && qryLen <= 2)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25))
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
                                    // indel-ins
                                    else if (refLen <= 2 && qryLen >= 3)
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                                    else
                                    {
                                        if ((abs(refStart-trueRefStart)<=1) &&
                                            (abs((refStart+refLen)-(trueRefStart+trueRefLen))<=10) &&
                                            ((abs(refLen-trueRefLen)/(float)trueRefLen)<=0.25) && 
                                            ((abs(qryLen-trueQryLen)/(float)trueQryLen)<=0.25))
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
                        if (genotypeTrueNum > 0)  // ����0�����ҵ���
                        {
                            // ���call��vcf-vector
                            call_length_list.push_back(trueSvLen);

                            true_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                        "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                            recall_length_list.push_back(trueSvLen);

                            // ������ȷ
                            if (genotypeTrueNum >= gtVec.size())  // ���ڵ���gtVec.size()���������ȷ
                            {
                                genotype_length_list.push_back(trueSvLen);
                                true_Gt_txt += "recall_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                            "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";
                                // ��ӳ��ȵ�rocTrueMap��
                                rocRecallMap[rocNum].push_back(svLength);
                            }
                            // ���ʹ���
                            else
                            {
                                genotypeMisTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n" + 
                                                "true_length\t" + to_string(trueSvLen) + "\t" + truevcfInfo + "\n";

                                // ���ж�����recall������roc�������������ֻ��genotype����roc
                                if (rocType == "recall")
                                {
                                    // ��ӳ��ȵ�rocTrueMap��
                                    rocRecallMap[rocNum].push_back(svLength);
                                }
                            }

                            // ɾ���Ѿ����ù��ı���
                            if (trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].find(refStart) != trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].end())
                            {
                                trueVcfStructure.chrStartLenInfoGtVecMap[chromosome].erase(refStart);
                            }
                        }
                        else  // ����ϱ�ѭ��û�ҵ��漯��ı��죬��������Լ��ҵ��ĳ��������
                        {
                            call_length_list.push_back(svLength);
                            misCallTxt += "call_length\t" + to_string(svLength) + "\t" + informations + "\n";
                        }
                    }
                    // ����ַ���
                    informations.clear();
                    string().swap(informations);

                    // ������
                    if (true_txt.size() > 10000000 || 
                        true_Gt_txt.size() > 10000000 || 
                        genotypeMisTxt.size() > 10000000 || 
                        misCallTxt.size() > 10000000) // ÿ10Mbд��һ��
                    {
                        gzwrite(trueFile, true_txt.c_str(), true_txt.length());
                        gzwrite(trueGtFile, true_Gt_txt.c_str(), true_Gt_txt.length());
                        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // ���ʹ���
                        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // �����漯�еı���

                        // ����ַ���
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
        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        // ������
        gzwrite(trueFile, true_txt.c_str(), true_txt.length());

        gzwrite(trueGtFile, true_Gt_txt.c_str(), true_Gt_txt.length());

        gzwrite(genotypeMisFile, genotypeMisTxt.c_str(), genotypeMisTxt.length()); // ���ʹ���

        gzwrite(misCallFile, misCallTxt.c_str(), misCallTxt.length()); // �����漯�еı���

        // ����û���ҵ����죬������
        saveFailCall(
            trueVcfStructure.chrStartLenInfoGtVecMap, 
            "failcall.vcf.gz"
        );

        // �ر��ļ�
        gzclose(trueFile);
        gzclose(trueGtFile);
        gzclose(genotypeMisFile);
        gzclose(misCallFile);

        int all_num = trueVcfStructure.allLengthList.size();
        int call_num = call_length_list.size();
        int recall_num = recall_length_list.size();

        int genotype_recall_number = genotype_length_list.size();

        // ������
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

        // ���Ƚ�����
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

        // �ͷ��ڴ�
        outFile.close();

        // ����ROC
        roc_calculate(rocCallMap, rocRecallMap, trueVcfStructure.allLengthList);

        return 0;
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
	 * ��ȡ����ĳ�����Ϣ
	 *
     * @param refLen            ref����
	 * @param qryLenVec         qry�����б�
     * @param gtVec             �������б�
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
	 * ͳ�Ʊ��쳤����Ϣ
	 *
     * @param sv_length         ���ֵ�����
	 * @param length_list       �����б�
     * 
     * 
     * @return vector<int>      ÿ������ĳ���
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
	 * ����recall��failCall�Ľ��
	 *
	 * @param chrStartLenInfoGtVecMap     �漯����ʣ�µ�vcf��Ϣ
     * @param outFileName                 ����ļ���
     * 
     * 
     * @return int             0
	**/
    int saveFailCall(
        map<string, map<int, tuple<int, vector<int>, string, int, vector<int> > > > & chrStartLenInfoGtVecMap, 
        const string & outFileName
    )
    {
        // û���ҵ��ı���
        string failCallTxt;
        // ����ļ���
        gzFile gzfp = gzopen(outFileName.c_str(), "wb");

        // �����ֵ䣬��û���ҵ���vcf���б���
        for (auto it1 : chrStartLenInfoGtVecMap)  // unordered_map<chr, map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> > >
        {
            for (auto it2 : it1.second)  // map<refStart, tuple<refLen, qryLenVec, vcfInfo, svLen> >
            {
                failCallTxt += get<2>(it2.second) + "\n";

                // ������
                if (failCallTxt.size() > 10000000) // ÿ10Mbд��һ��
                {
                    gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

                    // ����ַ���
                    failCallTxt.clear();
                    string().swap(failCallTxt);
                }
            }
        }
        // ������
        gzwrite(gzfp, failCallTxt.c_str(), failCallTxt.length());

        // �ͷ��ڴ棬�ر��ļ�
        gzclose(gzfp);

        return 0;
    }


    /**
	 * ����recall��failCall�Ľ��
	 *
	 * @param rocCallMap
     * @param rocRecallMap
     * @param all_length_list     ���б��쳤��
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
        // ������
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
        string outRoc = rocType + "." + rocKey + "\tTrue positives\tFalse positives\tTrue negatives\tFalse negatives\trecall\tprecision\tF-score\n";
        int allTruePositives = all_length_list.size(); // ������ȷ��λ������
        int truePositives = 0; // ������
        int falsePositives = 0; // ������
        int trueNegatives = 0; // ������
        int falseNegatives = 0; // ������

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
        allTruePositives = 0; // ������ȷ��λ������
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] == 0)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // ������
        falsePositives = 0; // ������
        trueNegatives = 0; // ������
        falseNegatives = 0; // ������

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all�е�snp����
            // ������score�µĳ���vector
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
                // ������score�µĳ���vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] == 0)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // ���������0��������score
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
        allTruePositives = 0; // ������ȷ��λ������
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] > -50 && all_length_list[i] < 50 && all_length_list[i] != 0)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // ������
        falsePositives = 0; // ������
        trueNegatives = 0; // ������
        falseNegatives = 0; // ������

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all�е�snp����
            // ������score�µĳ���vector
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
                // ������score�µĳ���vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] > -50 && findIter1->second[i] < 50 && findIter1->second[i] != 0)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // ���������0��������score
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
        allTruePositives = 0; // ������ȷ��λ������
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] <= -50)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // ������
        falsePositives = 0; // ������
        trueNegatives = 0; // ������
        falseNegatives = 0; // ������

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all�е�snp����
            // ������score�µĳ���vector
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
                // ������score�µĳ���vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] <= -50)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // ���������0��������score
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
        allTruePositives = 0; // ������ȷ��λ������
        for (size_t i = 0; i < all_length_list.size(); i++)
        {
            if (all_length_list[i] >= 50)
            {
                allTruePositives++;
            }
        }
        truePositives = 0; // ������
        falsePositives = 0; // ������
        trueNegatives = 0; // ������
        falseNegatives = 0; // ������

        for (auto iter1 = rocCallMap.rbegin(); iter1 != rocCallMap.rend(); iter1++)
        {
            int allPositiveTmp = 0;
            // all�е�snp����
            // ������score�µĳ���vector
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
                // ������score�µĳ���vector
                for (size_t i = 0; i < findIter1->second.size(); i++)
                {
                    if (findIter1->second[i] >= 50)
                    {
                        truePositiveTmp++;
                    }
                }
            }
            // ���������0��������score
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