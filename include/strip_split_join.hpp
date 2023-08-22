#ifndef STRIP_SPLIT_JOIN_HPP
#define STRIP_SPLIT_JOIN_HPP
#include <string>
#include <vector>

using namespace std;


/*
 Remove special characters from both ends of a string
 str -> string
 ch -> special character
*/
string strip(const string & str, char ch=' ');


/*
 Split the string
 str -> String
 delim -> Split characters
*/
vector<string> split(const string & str, const string & delim);


/*
 String join merging
 val -> The container to be joined
 delim -> Characters added between container elements
*/
string join(vector<int> & val, string delim);


/*
 String join merging
 val -> The container to be joined
 delim -> Characters added between container elements
*/
string join(vector<long int> & val, string delim);


/*
 String join merging
 val -> The container to be joined
 delim -> Characters added between container elements
*/
string join(vector<float> & val, string delim);


/*
 String join merging
 val -> The container to be joined
 delim -> Characters added between container elements
*/
string join(vector<string> & val, string delim);

#endif