#include "../include/strip_split_join.hpp"

using namespace std;

// Remove the character ch from both ends of str.
string strip(const string & str, char ch)
{
    const size_t begin = str.find_first_not_of(ch);
    if (begin == string::npos)
        return string();
    const size_t end = str.find_last_not_of(ch);
    return str.substr(begin, end + 1 - begin);
}

// Split str at every delim; empty fields are kept.
vector<string> split(const string & str, const string & delim)
{
    vector<string> res;
    split(str, delim, res);
    return res;
}

// Split str at every delim into fields, reusing the strings already in fields; empty fields are kept.
void split(const string & str, const string & delim, vector<string> & fields)
{
    size_t count = 0;
    if (!str.empty()) {
        if (delim.empty()) {
            fields.resize(1);
            fields[0] = str;
            return;
        }

        // ------ Count the fields so the vector grows at most once ------
        size_t fieldNum = 1;
        for (size_t pos = str.find(delim); pos != string::npos; pos = str.find(delim, pos + delim.size())) {
            ++fieldNum;
        }
        fields.reserve(fieldNum);

        size_t start = 0;
        while (true) {
            const size_t pos = str.find(delim, start);
            const size_t end = (pos == string::npos) ? str.size() : pos;
            if (count < fields.size()) {
                fields[count].assign(str, start, end - start);
            } else {
                fields.emplace_back(str, start, end - start);
            }
            ++count;
            if (pos == string::npos) break;
            start = pos + delim.size();
        }
    }
    fields.resize(count);
}

// Join integers with delim.
string join(const vector<int> & val, const string & delim)
{
    std::string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto& iter : val)
    {
        str += to_string(iter);

        if (index != vecSize-1)
        {
            str += delim;
        }
        index++;
    }
    return str;
}

// Join long integers with delim.
string join(const vector<long int> & val, const string & delim)
{
    std::string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto& iter : val)
    {
        str += to_string(iter);

        if (index != vecSize-1)
        {
            str += delim;
        }
        index++;
    }
    return str;
}

// Join floats (std::to_string) with delim.
string join(const vector<float> & val, const string & delim)
{
    std::string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto& iter : val)
    {
        str += to_string(iter);

        if (index != vecSize-1)
        {
            str += delim;
        }
        index++;
    }
    return str;
}

// Join strings with delim.
string join(const vector<string> & val, const string & delim)
{
    std::string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto& iter : val)
    {
        str += iter;

        if (index != vecSize-1)
        {
            str += delim;
        }
        index++;
    }
    return str;
}
