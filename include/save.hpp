#ifndef SAVE_HPP
#define SAVE_HPP

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

/**
 * @brief save result
 * 
 * @param outputFileName
 * 
 * @return
**/
class SAVE
{
private:
    // output stream
    string outputFileName_;

    // file handle
    ofstream fpO;
    gzFile gzfpO;
public:
    SAVE() {}
    SAVE(string aliFileName) {
        outputFileName_ = aliFileName;

        if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
        {
            // open file
            gzfpO = gzopen(outputFileName_.c_str(), "wb");
            if(!gzfpO)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << outputFileName_ << "': No such file or directory or possibly reached the maximum open file limit. You can set 'ulimit -n' to a larger value to continue." << endl;
                exit(1);
            }
        }
        else if (outputFileName_.size() > 0)
        {
            // open file
            fpO.open(outputFileName_, ios::out);
            if(!fpO)
            {
                cerr << "[" << __func__ << "::" << getTime() << "] " 
                    << "'" << outputFileName_ << "': No such file or directory or possibly reached the maximum open file limit. You can set 'ulimit -n' to a larger value to continue." << endl;
                exit(1);
            }
        }
    }
    ~SAVE() {
        if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
        {
            // close file
            gzclose(gzfpO);
        }
        else if (outputFileName_.size() > 0)
        {
            // close file
            fpO.close();
        }
    }
    
    int save(string & outTxt);
};

#endif