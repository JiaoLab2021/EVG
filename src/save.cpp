#include "../include/save.hpp"

using namespace std;

/**
 * @brief ��ʼ��
 * 
 * @param outputFileName_
 * 
 * @return
**/
void SAVE::SAVE::init(
    const string & outputFileName_
)
{
    outputFileName = outputFileName_;
}

/**
 * @brief ���ļ�
 * 
 * @return 0
**/
int SAVE::SAVE::open(

)
{
    if (outputFileName.find(".gz") != string::npos || outputFileName.find(".GZ") != string::npos)
    {
        // ����ļ���
        gzfpO = gzopen(outputFileName.c_str(), "wb");
        if(!gzfpO)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << outputFileName << "': No such file or directory." << endl;
            exit(1);
        }
    }
    else if (outputFileName.size() > 0)
    {
        // ����ļ���
        fpO.open(outputFileName, ios::out);
        if(!fpO)
        {
            cerr << "[" << __func__ << "::" << getTime() << "] " 
                << "'" << outputFileName << "': No such file or directory." << endl;
            exit(1);
        }
    }

    return 0;
}


/**
 * @brief �洢
 * 
 * @return int
**/
int SAVE::SAVE::save(
    string & outTxt_
)
{
    outTxt_ = strip(outTxt_, '\n');  // ȥ�����з�

    if (outputFileName.find(".gz") != string::npos || outputFileName.find(".GZ") != string::npos)
    {
        gzwrite(gzfpO, outTxt_.c_str(), outTxt_.length());
    }
    else if (outputFileName.size() > 0)
    {
        fpO << outTxt_ << endl;
    }
    else
    {
        cout << outTxt_ << endl;
    }
    
    return 0;
}

/**
 * @brief �ر��ļ�
 * 
 * @return 0
**/
int SAVE::SAVE::close(

)
{
    if (outputFileName.find(".gz") != string::npos || outputFileName.find(".GZ") != string::npos)
    {
        // �ر��ļ�
        gzclose(gzfpO);
    }
    else if (outputFileName.size() > 0)
    {
        // �ر��ļ�
        fpO.close();
    }

    return 0;
}