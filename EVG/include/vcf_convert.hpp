#ifndef vcf_convert_hpp
#define vcf_convert_hpp

#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <algorithm>
#include "zlib.h"
#include "strip_split_join.hpp"
#include "kseq.h"
#include "get_time.hpp"
#include "check_vcf_sort.hpp"
#include "save.hpp"
#include "vcf_open.hpp"

using namespace std;

namespace CONVERT{
    // kseq.h ���ļ�
    KSEQ_INIT(gzFile, gzread)

    struct refIndexStruct
    {
        map<string, string> sequenceMap;  // map<chromosome, sequence>
        string chrLenTxt;  // �洢contig���ȣ����浽vcfͷ
    };

    /*
        ��fasta�ļ�
        inputFileName -> refgenome
        refIndexS -> ����contig���Ⱥ�������Ϣ
    */
    void build_reference_index(string inputFileName, refIndexStruct & refIndexS)
    {
        // Ⱦɫ�����������graphtyper��Ҫ
        ofstream outFile;
        outFile.open("CHROMOSOME.NAME", ios::out);
        if(!outFile.is_open())
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                 << "'CHROMOSOME.NAME': No such file or directory." 
                 << endl;
            exit(1);
        }

        // �����ļ���
        gzFile gzfp = gzopen(inputFileName.c_str(), "rb");

        // ���ļ�
        if(!gzfp)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                 << "'"
                 << inputFileName 
                 << "': No such file or directory." 
                 << endl;
            exit(1);
        }
        else
        {
            kseq_t *ks;
            ks = kseq_init(gzfp);
        
            while( kseq_read(ks) >= 0 )
            {
                string chromosome = ks->name.s;
                long long int chrLen = ks->seq.l;
                string sequence = ks->seq.s;

                refIndexS.chrLenTxt += "##contig=<ID=" + chromosome + ",length=" + to_string(chrLen) + ">\n";
                refIndexS.sequenceMap[chromosome] = sequence;

                // ���Ⱦɫ������
                outFile << chromosome + "\n";
            }

            // �ͷ��ڴ棬�ر��ļ�
            kseq_destroy(ks);
            gzclose(gzfp);
        }

        // �ͷ��ڴ�
        outFile.close();
    }


    /*
        ��vcf�ļ������滻
        vcfFilename -> ��Ҫת����vcf�ļ�
        readLen -> �����ļ���read����
        refIndexS -> contig���Ⱥ�������Ϣ
        outFilename -> ����ļ���
        a. ��ͷ��Ⱦɫ�峤��
        b. ��һ������Ҫ����read length
        c. ���б�����reference��Ӧ
        d. ��vcf�ĵڰ��� 'END=' �����滻
        e. refSeq�ĵ�һ�����Ҫ��qrySeq�ĵ�һ�����һ��
        f. ���ref��qry�����Ƿ�һ����һ����λ������
        g. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
        h. ����滻���qry��û����ͬ�ģ��еĻ�������λ�� --> e.
        i. ����Ƿ���λ���ظ��ı���
        j. ���������е�.תΪ.|.
        k. ���������е�/תΪ|
        l. ֻ�����������еĶ��������
        m. ���GT�ǲ��Ǳ�qry�����л���
    */
    /**
     * @brief ��vcfת��Ϊ graph genome tools ��Ҫ�ĸ�ʽ
     * 
     * @param vcfFileName     ����vcf�ļ�
     * @param readLen         ����
     * @param MAF            �ε�λ����Ƶ��
     * @param MISSRATE       ȱʧ��
     * @param refIndexS       �ο������鳤�Ⱥ�������Ϣ
     * @param outFilename     ����ļ���
     * 
     * @return 0
    **/
    void vcf_convert(
        const string & vcfFilename, 
        const int & readLen, 
        const double & MAF, 
        const double & MISSRATE, 
        const refIndexStruct & refIndexS, 
        const string & outFilename
    )
    {
        // ����ļ��Ƿ�����
        check_vcf_sort(vcfFilename);

        // �������ϵ�������Ϣ
        map<string, string> seqMap = refIndexS.sequenceMap;

        // ����ת����Ľ��
        string outTxt;


        // ����ļ���
        SAVE::SAVE SAVECLASS;
        SAVECLASS.init(
            outFilename
        );
        SAVECLASS.open();


        // ��¼��һ�������start��Ⱦɫ���
        string preChromosome;
        long long int preRefStart = 0;


        // �����ļ���
        // �洢vcf��Ϣ
        VCFOPEN::VCFINFOSTRUCT INFOSTRUCTTMP;
        VCFOPEN::VCFOPEN VCFOPENCLASS;
        VCFOPENCLASS.init(
            vcfFilename
        );
        // ��vcf�ļ�
        VCFOPENCLASS.open();


        // ���û�б����꣬����
        while (VCFOPENCLASS.read(INFOSTRUCTTMP))
        {
            // ע����
            if (INFOSTRUCTTMP.INFO.find("#") != string::npos)
            {
                // ���vcf�к���Ⱦɫ�峤����Ϣ������
                if (INFOSTRUCTTMP.INFO.find("#CHROM") != string::npos)
                {
                    // ����Ⱦɫ�峤����Ϣ
                    outTxt += refIndexS.chrLenTxt;

                    // �����ͷ
                    outTxt += INFOSTRUCTTMP.INFO + "\n";
                }
                else if (INFOSTRUCTTMP.INFO.find(",length") == string::npos)
                {
                    outTxt += INFOSTRUCTTMP.INFO + "\n";
                }

                continue;
            }

            // ��ע����
            if (INFOSTRUCTTMP.INFOVec.size() < 9) // �ȼ���ļ��Բ��ԣ���������������
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error '"
                        << vcfFilename 
                        << "' -> The number of vcf columns is less than 9." 
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;
                exit(1);
            }



            /* ************************ filter SNPs by maf and missing rate ************************ */
            // �ж��Ƿ���Ҫ����
            if (MAF > 0 && MISSRATE < 1)
            {
                // ��ȡ��������
                INFOSTRUCTTMP.TYPE = VCFOPENCLASS.get_TYPE(
                    INFOSTRUCTTMP.LEN, 
                    INFOSTRUCTTMP.ALTVec
                );

                if (INFOSTRUCTTMP.TYPE == "SNP")  // SNPʱ���ж��Ƿ����
                {
                    double MAFTMP;  // ��С��λ����Ƶ��
                    double MISSRATETMP;  // ȱʧ��

                    // ��ȡ���еĻ�����   map<idx, vector<gtString>>
                    map<int, vector<string> > GTVecMapTmp = VCFOPENCLASS.get_gt(
                        INFOSTRUCTTMP.INFOVec
                    );

                    // ���ֻ��һ�������ͣ�������
                    if (GTVecMapTmp.size() <= 1)
                    {
                        continue;
                    }
                    
                    // ������С��λ����Ƶ�ʺ�ȱʧ��
                    tie(MAFTMP, MISSRATETMP) = VCFOPENCLASS.calculate(
                        GTVecMapTmp, 
                        INFOSTRUCTTMP.INFOVec.size() - 9
                    );

                    // ûͨ����ֵ��ֱ����һ��ѭ��
                    if (MAFTMP < MAF || MISSRATETMP > MISSRATE)
                    {
                        continue;
                    }
                }
            }



            // �������ϵ�ref������Ϣ
            // �ȿ��ύ�Ļ���������û�ж�Ӧ��Ⱦɫ����Ϣ��û�еĻ��˳�����
            if (seqMap.find(INFOSTRUCTTMP.CHROM) == seqMap.end())
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error: "
                        << INFOSTRUCTTMP.CHROM 
                        << " -> Chromosome does not exist in refgenome." 
                        << endl;
                exit(1);
            }

            // ���qrySeq��û��<֮�����Ϣ���еĻ���λ��ֱ��д�� (<INS>;<DUP>)
            if (INFOSTRUCTTMP.INFOVec[4].find("<") != string::npos)
            {
                outTxt += INFOSTRUCTTMP.INFO + "\n";

                continue;
            }
            

            // 1. ����һ�������Ƿ����read length��������ڵĻ��ٽ���ת�������������ñ���
            if (readLen > INFOSTRUCTTMP.POS) // paragraph��Ҫ��
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: start of variation is less than read length, skip this site -> "
                        << INFOSTRUCTTMP.CHROM << "\t" 
                        << INFOSTRUCTTMP.POS << endl;

                continue;
            }
            

            // 2. vcf��Ӧ��reference�ϵ�����
            if (seqMap[INFOSTRUCTTMP.CHROM].size() < (INFOSTRUCTTMP.END)) // ���Ⱦɫ�峤���Ƿ���ȷ
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Chromosome length error: " 
                        << INFOSTRUCTTMP.CHROM << ".size()=" << seqMap[INFOSTRUCTTMP.CHROM].size() << endl;
                exit(1);
            }
            
            string trueRefSeq = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN);
            if (trueRefSeq != INFOSTRUCTTMP.REF) // ����Ͳο��������ϵ����в�һ�������滻�ɲο������������
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: sequence difference between refgenome and vcf, replace by refgenome sequence -> "
                        << INFOSTRUCTTMP.CHROM << " " 
                        << INFOSTRUCTTMP.POS << " " 
                        << INFOSTRUCTTMP.REF << "(" << trueRefSeq << ")\n";

                INFOSTRUCTTMP.INFOVec[3] = trueRefSeq;
                INFOSTRUCTTMP.REF = trueRefSeq;
            }


            // 3. ��vcf�ĵڰ��� 'END=' �����滻
            // paragrpah (raise Exception("{}:{} error in adding ref support.".format(start, end)))
            // ������ʽ��END���������滻
            std::regex endReg("END=\\d+");
            INFOSTRUCTTMP.INFOVec[7] = regex_replace(INFOSTRUCTTMP.INFOVec[7], endReg, "END=" + to_string(INFOSTRUCTTMP.END));


            // 4. Different padding base for REF and ALT. (paragraphҪ��SV��ALT��һ����ĸҪ��REFһ��)
            for (int i = 0; i < INFOSTRUCTTMP.ALTVec.size(); i++)
            {
                string qrySeq = INFOSTRUCTTMP.ALTVec[i];

                // ��һ�������һ����ʱ�򣬶�qrySeq�����м���refSrq��ǰһ�����
                if (qrySeq[0] != INFOSTRUCTTMP.REF[0] && (qrySeq.length() > 1 || INFOSTRUCTTMP.REF.length() > 1))
                {
                    // ������ǰŲ1
                    INFOSTRUCTTMP.POS = INFOSTRUCTTMP.POS - 1;
                    INFOSTRUCTTMP.INFOVec[1] = to_string(INFOSTRUCTTMP.POS);
                    INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.LEN + 1;
                    INFOSTRUCTTMP.REF = seqMap[INFOSTRUCTTMP.CHROM].substr(INFOSTRUCTTMP.POS - 1, INFOSTRUCTTMP.LEN); // ������ȡ������Ϣ
                    INFOSTRUCTTMP.INFOVec[3] = INFOSTRUCTTMP.REF; // ��Vector���и�ֵ

                    // �� 'INFOSTRUCTTMP.ALTVec' ���������ж����� 'refSeq[0]'
                    for (size_t j = 0; j < INFOSTRUCTTMP.ALTVec.size(); j++)
                    {
                        INFOSTRUCTTMP.ALTVec[j] = INFOSTRUCTTMP.REF[0] + INFOSTRUCTTMP.ALTVec[j];
                    }

                    qrySeq = INFOSTRUCTTMP.ALTVec[i]; // qrySeq ���¸�ֵ
                }

                // 5. ���refSeq��qrySeq�Ƿ�һ����һ��������
                if (qrySeq == INFOSTRUCTTMP.REF)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: sequence same in REF and ALT, skip this site -> "
                            << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                            << endl;

                    continue;
                }

                // 6. ���refSeq��qrySeq���Ƿ���atgcnATGCN����ַ������еĻ�������λ��
                smatch results;
                std::regex atgcReg("[^ATGCNatgcn]");
                if (regex_search(qrySeq, results, atgcReg) != 0)
                {
                    cerr << "[" << __func__ << "::" << getTime() << "] "
                            << "Warning: sequence contains non-ATGCNatgcn characters, skip this site -> "
                            << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                            << endl;

                    continue;
                }
            }
            // ��qry���н����滻
            INFOSTRUCTTMP.INFOVec[4] = join(INFOSTRUCTTMP.ALTVec, ",");


            // 7. ����滻���qry��û����ͬ�ģ��еĻ�������λ��
            // bayestyperҪ�� (A    ATG,TGT,ATG)
            // Assertion `count(alt_alleles.begin() + i + 1, alt_alleles.end(), alt_alleles.at(i)) == 0' failed.
            int qrySeqNum = 0;
            for (auto it : INFOSTRUCTTMP.ALTVec)
            {
                int qrySeqNumTmp = count(INFOSTRUCTTMP.ALTVec.begin(), INFOSTRUCTTMP.ALTVec.end(), it);
                qrySeqNum = max(qrySeqNumTmp, qrySeqNum);
            }
            if (qrySeqNum > 1)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: allelic repeat, skip -> "
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;

                continue;
            }


            // 8. �������Ƿ��ظ����ظ��Ļ�����
            // ������µ�Ⱦɫ�壬����ʼλ�ù���
            if (INFOSTRUCTTMP.CHROM != preChromosome)
            {
                preChromosome = INFOSTRUCTTMP.CHROM;
                preRefStart = 0;
            }
            if (INFOSTRUCTTMP.POS == preRefStart)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: multiple variants observed on position " 
                        << INFOSTRUCTTMP.POS 
                        << " on " 
                        << INFOSTRUCTTMP.CHROM
                        << ", skip this site"
                        << endl;

                continue;
            }


            // 9. �Ի������н��м�顣ֻ���ж�������졢'/'�滻Ϊ'|'��'.'�滻Ϊ'.|.'��PanGenie���
            // ��FORMAT�ֶ���gt��λ��
            vector<string> formatVec = split(strip(INFOSTRUCTTMP.INFOVec[8], '\n'), ":");
            vector<string>::iterator gtItera = find(formatVec.begin(), formatVec.end(), "GT");
            int maxGT = 0; // ��¼����GT�����Ƿ��qry�����л��࣬��Ļ�������λ��
            // ֻҪGT�ֶ�
            INFOSTRUCTTMP.INFOVec[8] = "GT";

            int gtIndex = 0;
            if (gtItera != formatVec.end()) // FORMAT����GT
            {
                // GT��index
                gtIndex = distance(formatVec.begin(), gtItera);

                // ��GT�н���ѭ��
                for (int i = 9; i < INFOSTRUCTTMP.INFOVec.size(); i++)
                {
                    // �ҵ�������
                    string gt = strip(split(INFOSTRUCTTMP.INFOVec[i], ":")[gtIndex], '\n');

                    // ���ǲ���'.'���ǵĻ��滻Ϊ'.|.'
                    if (gt == ".")
                    {
                        gt = ".|.";
                    }

                    // ��'/'�滻Ϊ'|'
                    if (gt.find("/") != string::npos)
                    {
                        std::regex reg("/");
                        gt = regex_replace(gt, reg, "|");
                    }
                    INFOSTRUCTTMP.INFOVec[i] = gt;

                    // ѭ��������GT
                    if (gt.find("/") != string::npos)
                    {
                        for (auto it1 : split(gt, "/"))
                        {
                            try
                            {
                                if (stoi(it1) > maxGT)
                                {
                                    maxGT = stoi(it1);
                                }
                            }
                            catch(const std::invalid_argument& e)
                            {
                                continue;
                            }
                        }
                    }
                    else if (gt.find("|") != string::npos)
                    {
                        for (auto it1 : split(gt, "|"))
                        {
                            try
                            {
                                if (stoi(it1) > maxGT)
                                {
                                    maxGT = stoi(it1);
                                }
                            }
                            catch(const std::invalid_argument& e)
                            {
                                continue;
                            }
                        }
                    }
                    else
                    {
                        cerr << "[" << __func__ << "::" << getTime() << "] "
                                << "Error: invalid delimiter: " 
                                << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                                << endl;
                        exit(1);
                    }
                }
            }
            else // û�еĻ��˳�����
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Error: GT not in FORMAT column -> " 
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;
                exit(1);
            }


            // 10. ���GT�ǲ��Ǳ�qry�����л���
            // pangenieҪ��VariantReader::VariantReader: invalid genotype in VCF.
            if (INFOSTRUCTTMP.ALTVec.size() < maxGT)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] "
                        << "Warning: there are more GTs than qry sequences: " 
                        << " (" << strip(INFOSTRUCTTMP.INFO, '\n') << ")"
                        << endl;

                continue;
            }
            

            // ��preRefStart���¸�ֵ
            preRefStart = INFOSTRUCTTMP.POS;

            // ���information�Ƿ���գ���յĻ�����
            if (INFOSTRUCTTMP.INFO.empty())
            {
                continue;
            }

            // ���滻����ַ����ӵ�outTxt��
            outTxt += join(INFOSTRUCTTMP.INFOVec, "\t") + "\n";

            if (outTxt.size() > 10000000) // ÿ10Mbд��һ�Σ����ٴ���IO
            {
                SAVECLASS.save(
                    outTxt
                );

                // ����ַ���
                outTxt.clear();
                string().swap(outTxt);
            }
        }


        if (outTxt.size() >= 0)  // ���дһ��
        {
            SAVECLASS.save(
                outTxt
            );

            // ���
            outTxt.clear();
            string().swap(outTxt);
        }

        // �ر��ļ�
        VCFOPENCLASS.close();
        SAVECLASS.close();
    }

}  // namespace CONVERT

#endif