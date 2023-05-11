#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <getopt.h>
#include "get_time.hpp"
#include "strip_split_join.hpp"

using namespace std;

namespace convert
{
    int convert(string * fastaFilename, string * outputFilename)
{
    ifstream fastaFile;
    fastaFile.open(*fastaFilename, ios::in);

    string informations;
    string outTxt;

    if(!fastaFile.is_open())
    {
        cerr << "[" << __func__ << "::" << getTime() << "] '" 
                 << *fastaFilename << "': No such file or directory.\n";
             exit(1);
    }
    else
    {
        while(getline(fastaFile,informations))
        {
            if(informations.empty())
            {
                continue;
            }
            else
            {   
                string::size_type idx = informations.find(">");
                if ( idx != string::npos )
                {
                    outTxt += "\n" + informations + "\n";
                }
                else
                {
                    outTxt += informations;
                }
            }
        }
    }

    // 关闭文件
    fastaFile.close();

    ofstream outputFile;
    outputFile.open(*outputFilename);

    outputFile << strip(outTxt, '\n') << "\n";

    // 关闭文件
    outputFile.close();

    return 0;
}
}