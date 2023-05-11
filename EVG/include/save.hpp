#ifndef save_hpp
#define save_hpp

#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <map>
#include <tuple>
#include <iomanip>
#include "zlib.h"
#include <unordered_map>
#include "strip_split_join.hpp"
#include "get_time.hpp"

using namespace std;

namespace SAVE
{
    /**
     * @brief 打开vcf文件
     * 
     * @param outputFileName_
     * 
     * @return
    **/
    class SAVE
    {
    private:
        // vcf文件
        string outputFileName;

        // 输出文件流
        ofstream fpO;
        gzFile gzfpO;
    public:
        void init(
            const string & outputFileName_
        );

        int open();
        int save(
            string & outTxt_
        );
        int close();
    };

    /**
     * @brief 初始化
     * 
     * @param outputFileName_
     * 
     * @return
    **/
    void SAVE::init(
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
   int SAVE::open(

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
    int SAVE::save(
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
    int SAVE::close(

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

} // namespace SAVE

#endif