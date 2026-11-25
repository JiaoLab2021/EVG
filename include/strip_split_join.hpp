#pragma once

#include <string>
#include <vector>

using namespace std;

// Remove the character ch from both ends of str.
string strip(const string & str, char ch=' ');

// Split str at every delim; empty fields are kept. The second form fills fields, reusing the strings already in it.
vector<string> split(const string & str, const string & delim);
void split(const string & str, const string & delim, vector<string> & fields);

// Join numbers (std::to_string) or strings with delim.
string join(const vector<int> & val, const string & delim);
string join(const vector<long int> & val, const string & delim);
string join(const vector<float> & val, const string & delim);
string join(const vector<string> & val, const string & delim);
