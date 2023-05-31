#ifndef STRIP_SPLIT_JOIN_HPP
#define STRIP_SPLIT_JOIN_HPP
#include <string>
#include <vector>

using namespace std;


/*
	去掉字符串两端的特殊字符
	str -> 字符串
	ch -> 特殊字符
*/
string strip(const string & str, char ch=' ');


/*
	对字符串进行拆分
	str -> 字符串
	delim -> 拆分的字符
*/
vector<string> split(const string & str, const string & delim);


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<int> & val, string delim);


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<long int> & val, string delim);


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<float> & val, string delim);


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<string> & val, string delim);

#endif