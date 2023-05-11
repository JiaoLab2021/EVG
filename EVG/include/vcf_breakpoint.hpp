#ifndef vcf_breakpoint_hpp
#define vcf_breakpoint_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

namespace BREAKPOINT
{
    int vcf_breakpoint(const string & vcfFileName, const string & prefix, const int & breakpointErrorSize)
    {
        // �����ļ���
        gzFile gzfpI = gzopen(vcfFileName.c_str(), "rb");
        if(!gzfpI)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " << "'" << vcfFileName << "': No such file or directory." << endl;
            exit(1);
        }

        // ��¼ת����Ľ��
        string outTxt;

        // ����ļ���
        string outFileName = prefix + ".vcf.gz";

        // ����ļ���
        gzFile outFile = gzopen(outFileName.c_str(), "wb");

        if(!outFile)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "'"
                 << outFileName
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }

        // ��vcf�ļ�����ת��
        string information = ""; // ��ʱ�ַ���
        char line[1024]; // һ��ֻ��1024�ֽڵ�����

        while(gzgets(gzfpI, line, 1024))
        {
            information += line;
            if (information.find("\n") != string::npos) // һ�н���
            {
                if(line != nullptr && line[0] == '\0')
                {
                    continue;
                }

                information = strip(information, '\n'); // ȥ�����з�

                if (information.find("#") != string::npos) // ����Ƿ���ע����
                {
                    string headLine = information + "\n";
                    gzwrite(outFile, headLine.c_str(), headLine.length());
                }
                else
                {
                    vector<string> informationVec = split(information, "\t");

                    // vcf�ĸ�����Ϣ
                    string refChr = informationVec[0];
                    long long int refStart = stol(informationVec[1]);

                    string refSeq = informationVec[3];
                    string qrySeqs = informationVec[4];
                    vector<string> qrySeqsVec = split(qrySeqs, ",");
                    long long int refLen = refSeq.length();
                    long long int refEnd = refStart + refLen - 1;
                    long long int svLen = qrySeqs.length() - refLen;


                    // ����INFO��ߵ�end��Ϣ
                    smatch endResult;
                    regex endReg("END=(\\d+)");

                    string::const_iterator iterStart = informationVec[7].begin();
                    string::const_iterator iterEnd = informationVec[7].end();
                    regex_search(iterStart, iterEnd, endResult, endReg);

                    string endResultTmp = endResult[1];
                    if (endResultTmp.empty())
                    {
                        if (refEnd == stol(endResult[1]))
                        {
                            continue;
                        }
                        else
                        {
                            cerr << "[" << __func__ << "::" << getTime() << "] " <<  information << endl;
                        }
                    }

                    // ������
                    refStart -= breakpointErrorSize;
                    informationVec[1] = to_string(refStart);
                    refEnd -= breakpointErrorSize;
                    informationVec[7] = regex_replace(informationVec[7], endReg, "END=" + to_string(refEnd));

                    // ��ӵ�����ַ�����
                    for (size_t i = 0; i < informationVec.size(); i++)
                    {
                        if (i > 0)
                        {
                            outTxt += "\t";
                        }
                        outTxt += informationVec[i];
                    }
                    outTxt += "\n";

                    // ������
                    if (outTxt.size() > 10000000) // ÿ10Mbд��һ��
                    {
                        gzwrite(outFile, outTxt.c_str(), outTxt.length());

                        // ����ַ���
                        outTxt.clear();
                        string().swap(outTxt);
                    }
                }

                // ����ַ���
                information.clear();
                string().swap(information);
            }
        }
        // �ر��ļ�
        gzclose(gzfpI);

        // ���д��һ��
        gzwrite(outFile, outTxt.c_str(), outTxt.length());

        // ����ַ���
        outTxt.clear();
        string().swap(outTxt);

        gzclose(outFile);

        return 0;
    }
} // namespace BREAKPOINT

#endif