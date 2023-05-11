#ifndef get_time_hpp
#define get_time_hpp
#include <string>
#include <time.h>

using namespace std;

string getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
    return tmp;
}

#endif