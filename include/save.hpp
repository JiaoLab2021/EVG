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

} // namespace SAVE

#endif