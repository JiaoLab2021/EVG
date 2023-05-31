#include "../include/save.hpp"

using namespace std;

/**
 * @brief 初始化
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
 * @brief 打开文件
 * 
 * @return 0
**/
int SAVE::SAVE::open(

)
{
    if (outputFileName.find(".gz") != string::npos || outputFileName.find(".GZ") != string::npos)
    {
        // 输出文件流
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
        // 输出文件流
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
 * @brief 存储
 * 
 * @return int
**/
int SAVE::SAVE::save(
    string & outTxt_
)
{
    outTxt_ = strip(outTxt_, '\n');  // 去除换行符

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
 * @brief 关闭文件
 * 
 * @return 0
**/
int SAVE::SAVE::close(

)
{
    if (outputFileName.find(".gz") != string::npos || outputFileName.find(".GZ") != string::npos)
    {
        // 关闭文件
        gzclose(gzfpO);
    }
    else if (outputFileName.size() > 0)
    {
        // 关闭文件
        fpO.close();
    }

    return 0;
}