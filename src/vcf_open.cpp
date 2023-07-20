#include "../include/vcf_open.hpp"

using namespace std;


/**
 * @brief ��ʼ��
 * 
 * @param vcfFileName   the output of vcf file
 * 
**/
VCFOPEN::VCFOPEN(
    const string & vcfFileName
)
{
    vcfFileName_ = vcfFileName;

    // �����ļ���
    gzfpI = gzopen(vcfFileName_.c_str(), "rb");
    if(!gzfpI)
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " 
            << "'" << vcfFileName_ << "': No such file or directory." << endl;
        exit(1);
    }
}

/**
 * @brief �ر��ļ�
 * 
 * @return 0
**/
VCFOPEN::~VCFOPEN()
{
    // �ر��ļ�
    gzclose(gzfpI);
}


// ��սṹ��
void VCFINFOSTRUCT::clear()
{
    // general information
    line.clear();
    lineVec.clear();
    ALTVec.clear();

    CHROM.clear();  //  [0]
    POS = 0;  //  [1]
    ID.clear();  // [2]
    REF.clear();  //  [3]
    ALT.clear();  // [4]
    QUAL = 0;  //  [5]
    FILTER.clear(); //  [6]
    INFO.clear();  // [7]
    FORMAT.clear();  // [8]

    LEN = 0;  // REF length
    END = 0;  // ALT length
}


/**
 * @brief traverse files
 * 
 * @param INFOSTRUCTTMP store the contents of the line
 * 
 * @return bool
**/
bool VCFOPEN::read(
    VCFINFOSTRUCT & INFOSTRUCTTMP
)
{
    INFOSTRUCTTMP.clear();

    string info = ""; // temporary string
    char line[1024]; // Read only 1024 bytes of data at a time

    if(gzgets(gzfpI, line, 1024))
    {
        info += line;

        // empty line
        memset(line, '\0', sizeof(line));

        // If there is no newline character, it means that the line is not over, continue to read
        while ((info.find("\n") == string::npos || info == "\n") && gzgets(gzfpI, line, 1024))
        {
            info += line;

            // empty line
            memset(line, '\0', sizeof(line));
        }
        
        // remove line breaks
        if (info.size() > 0)
        {
            info = strip(info, '\n');
        }
        else  // false
        {
            return false;
        }
        
        // assignment
        INFOSTRUCTTMP.line = info;
        INFOSTRUCTTMP.lineVec = split(info, "\t");  // segmentation

        // comment lines
        if (info.find("#") != string::npos)
        {
            return true;
        }
        

        // non-comment lines
        if (INFOSTRUCTTMP.lineVec.size() < 9) // Check the file first, if it is wrong, it will jump out of the code
        {
            cerr << "[" << __func__ << "::" << getTime() << "] "
                << "Error: '"
                << vcfFileName_ 
                << "' -> The number of vcf columns is less than 9.\n" 
                << INFOSTRUCTTMP.line << endl;
            exit(1);
        }

        INFOSTRUCTTMP.CHROM = INFOSTRUCTTMP.lineVec[0];  // CHROM
        INFOSTRUCTTMP.POS = stoul(INFOSTRUCTTMP.lineVec[1]);  // POS
        INFOSTRUCTTMP.ID = INFOSTRUCTTMP.lineVec[2];  // ID
        INFOSTRUCTTMP.REF = INFOSTRUCTTMP.lineVec[3];  // REF
        INFOSTRUCTTMP.ALT = INFOSTRUCTTMP.lineVec[4];  // ALT
        INFOSTRUCTTMP.ALTVec = split(INFOSTRUCTTMP.ALT, ",");  // 'ALT' �� 'vec'

        if (isdigit(INFOSTRUCTTMP.lineVec[5][0]))  // Assign QUAL if it is a number, otherwise 0.0
        {
            INFOSTRUCTTMP.QUAL = stod(INFOSTRUCTTMP.lineVec[5]);  // QUAL
        }
        else
        {
            INFOSTRUCTTMP.QUAL = 0.0;
        }

        INFOSTRUCTTMP.FILTER = INFOSTRUCTTMP.lineVec[6];  // FILTER
        INFOSTRUCTTMP.INFO = INFOSTRUCTTMP.lineVec[7];  // [7]
        INFOSTRUCTTMP.FORMAT = INFOSTRUCTTMP.lineVec[8];  // [8]

        INFOSTRUCTTMP.LEN = INFOSTRUCTTMP.REF.size();  // REF length
        INFOSTRUCTTMP.END = INFOSTRUCTTMP.POS + INFOSTRUCTTMP.LEN - 1;  // end position
    }
    else  // Return false directly after traversing
    {
        return false;
    }

    return true;
}


/**
 * Get the type of variant
 * 
 * @param refLen  the length of REF variant
 * @param ALTVec  vector<qrySeq>
 * 
 * @return string   TYPE: SNP, InDel, Deletion, Insertion, Inversion, Duplication, Other
**/
string VCFOPEN::get_TYPE(
    const uint32_t & refLen,
    const vector<string> & ALTVec
)
{
    uint32_t qryLen = 0;
    for (const auto& it1 : ALTVec)
    {
        uint32_t qryLenTmp = it1.size();
        qryLen = max(qryLen, qryLenTmp);  // Find the longest qrySeq
    }

    int32_t svLen = qryLen - refLen;
    double lengthRatio = qryLen / float(refLen);

    if (svLen == 0 && refLen == 1 && qryLen == 1)
    {
        return "SNP";
    }
    else if (svLen <= 49 && svLen >= -49 && refLen <= 49 && qryLen <= 49)
    {
        return "InDel";
    }
    else if (svLen >= -2 && svLen <= 2 && refLen > 49 && qryLen > 49 )
    {
        return "Inversion";
    }
    else if (lengthRatio >= 1.8 && lengthRatio <= 2.2 && refLen > 49 && qryLen > 49)
    {
        return "Duplication";
    }
    else if (svLen < 0)
    {
        return "Deletion";
    }
    else if (svLen > 0)
    {
        return "Insertion";
    }
    else
    {
        return "Other";
    }
    
    return "Other";
}


/**
 * Get a list of locus genotypes.
 *
 * @param lineTmpVec  lineTmpVec
 * 
 * 
 * @return GTVecMap   map<int, vector<string> >,  map<idx, vector<GTString> >
**/
map<int, vector<string> > VCFOPEN::get_gt(
    const vector<string> & lineTmpVec
)
{
    map<int, vector<string> > GTVecMap;  // ����Line���͵�map,  map<idx, vector<GTString> >

    // FORMAT������
    vector<string> formatVec = split(lineTmpVec[8], ":");  // FORMAT�ַ����

    int gtIndex = distance(formatVec.begin(), find(formatVec.begin(), formatVec.end(), "GT"));  // ��ȡGT������λ��

    if (gtIndex == formatVec.size())  // �ж�index�Ƿ���ڣ������ڵĻ����ػ����Ͷ�Ϊ0��
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Warning: no [GT] information in FORMAT -> " << lineTmpVec[0] << ":" << lineTmpVec[1] << endl;
        return GTVecMap;
    }
    else  // ������ڣ�����б���
    {
        string GTString;  // �洢�������ֶ�

        for (size_t i = 9; i < lineTmpVec.size(); i++)  // �ӵھ��п�ʼѭ��
        {
            GTString = split(lineTmpVec[i], ":")[gtIndex];  // gt�ֶ�

            // ���� '.' ��λ��ֱ������
            if (GTString.find(".") != string::npos)
            {
                continue;
            }
            
            string splitStr;  // gt�еķָ���
            if (GTString.find("/") != string::npos)  // �жϡ�/���ָ���
            {
                splitStr = "/";
            }
            else if (GTString.find("|") != string::npos)  // �жϡ�|��Ϊ�ָ���
            {
                splitStr = "|";
            }
            
            // ���֪���ָ���ٸ�ֵ
            if (splitStr.size() > 0)
            {
                GTVecMap[i - 9] = split(GTString, splitStr);  // ��ֵ
            }
        }
    }

    return GTVecMap;
}


/**
 * split gt
 *
 * @param gtTxt
 * 
 * 
 * @return vector<gt>
**/
vector<string> VCFOPEN::gt_split(const string & gtTxt)
{
    // ��ʱgt�б�
    vector<string> gtVecTmp;

    if (gtTxt.find("/") != string::npos)
    {
        gtVecTmp = split(gtTxt, "/");
    }
    else if (gtTxt.find("|") != string::npos)
    {
        gtVecTmp = split(gtTxt, "|");
    }
    else
    {
        cerr << "[" << __func__ << "::" << getTime() << "] " << "Error: GT is not separated by '/' or '|' -> " << gtTxt << endl;
        exit(1);
    }

    return gtVecTmp;
}


/**
 * Calculate MAF and MISSRATE.
 *
 * @param GTVecMap    map<int, vector<string> >,  map<idx, vector<GTString> >
 * @param sampleNum   total number of samples
 * 
 * 
 * @return tuple<double, double>   tuple<MAF, MISSRATE>
**/
tuple<double, double> VCFOPEN::calculate(
    const map<int, vector<string> > & GTVecMap, 
    uint32_t sampleNum
)
{
    double MAFTMP;  // ��С��λ����Ƶ��
    double MISSRATETMP;  // ȱʧ��

    /* ******************************** ���� MAF ******************************** */ 
    map<string, uint32_t> GTFreMapTmp;  // ��¼������Ƶ��  map<gt, fre>
    uint16_t ploidy = 1;  // ��¼����

    for (const auto& it1 : GTVecMap)  // map<idx, vector<GTString> >
    {
        for (const auto& it2 : it1.second)  // vector<GTString>
        {
            GTFreMapTmp[it2]++;  // �����Ͷ�Ӧ��Ƶ�ʼ�1
        }
        // �����1�ټ���
        if (ploidy == 1)
        {
            ploidy = it1.second.size();  //��¼����
        }
    }

    map<uint32_t, string> GTFreMapTmpTmp;  // mapת��, map<fre, gt>
    for (const auto& it1 : GTFreMapTmp)  // map<gt, fre>
    {
        GTFreMapTmpTmp[it1.second] = it1.first;
    }

    if (GTFreMapTmpTmp.size() > 1)
    {
        // �����ڶ��Ǵε�λ����
        auto iter1 = GTFreMapTmpTmp.end();
        iter1--; iter1--;
        // frequence/(ploidy*N)
        MAFTMP = iter1->first/(double)(ploidy*sampleNum);
    }
    
    /* ******************************** ���� MISSRATE ******************************** */ 
    // 1- number/allNumber
    MISSRATETMP = 1 - (GTVecMap.size()/(double)sampleNum);

    return make_tuple(MAFTMP, MISSRATETMP);
}