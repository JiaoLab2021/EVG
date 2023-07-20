#include "../include/save.hpp"

using namespace std;


/**
 * @brief save result to file
 * 
 * @return int
**/
int SAVE::save(
    string & outTxt
)
{
    outTxt = strip(outTxt, '\n');  // remove '\n'

    if (outTxt.empty())
    {
        return 0;
    }

    if (outputFileName_.find(".gz") != string::npos || outputFileName_.find(".GZ") != string::npos)
    {
        outTxt += "\n";
        gzwrite(gzfpO, outTxt.c_str(), outTxt.length());
    }
    else if (outputFileName_.size() > 0)
    {
        fpO << outTxt << endl;
    }
    else
    {
        cout << outTxt << endl;
    }
    
    return 0;
}