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
     * @brief ��vcf�ļ�
     * 
     * @param outputFileName_
     * 
     * @return
    **/
    class SAVE
    {
    private:
        // vcf�ļ�
        string outputFileName;

        // ����ļ���
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
     * @brief ��ʼ��
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
     * @brief ���ļ�
     * 
     * @return 0
    **/
   int SAVE::open(

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
    int SAVE::save(
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
    int SAVE::close(

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

} // namespace SAVE

#endif