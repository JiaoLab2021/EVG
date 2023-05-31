#include "../include/strip_split_join.hpp"

using namespace std;


/*
	去掉字符串两端的特殊字符
	str -> 字符串
	ch -> 特殊字符
*/
string strip(const string & str, char ch)
{
	int i = 0;
	while (str[i] == ch)  // 头部ch字符个数是 i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch)  // 结尾ch字符个数是 str.size() - 1 - j
		j--;		
	return str.substr(i, j+1-i);
}


/*
	对字符串进行拆分
	str -> 字符串
	delim -> 拆分的字符
*/
vector<string> split(const string & str, const string & delim)
{
	vector<string> res;  // 将分割后的子字符串存储在vector中
	if("" == str) return  res;  // 如果为空返回空vector
	
	string strs = str + delim;  // 扩展字符串以方便检索最后一个分隔出的字符串
	size_t pos;
	size_t size = strs.size();
 
	for (int i = 0; i < size; ++i) {
		pos = strs.find(delim, i); // pos为分隔符第一次出现的位置，从i到pos之前的字符串是分隔出来的字符串
		if( pos < size)  // 如果查找到，如果没有查找到分隔符，pos为string::npos
		{
			string s = strs.substr(i, pos - i);  // 从i开始长度为pos-i的子字符串
			res.push_back(s);  // 两个连续空格之间切割出的字符串为空字符串，这里没有判断s是否为空，所以最后的结果中有空字符的输出，
			i = pos + delim.size() - 1;
		}
	}
	return res;	
}


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<int> & val, string delim)
{
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val)
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


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<long int> & val, string delim)
{
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val)
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


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<float> & val, string delim)
{
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val)
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


/*
	字符串join合并
	val -> 需要join的容器
	delim -> 容器元素之间加入的字符
*/
string join(vector<string> & val, string delim)
{
    std::string str;
	int vecSize = val.size();
	int index = 0;
	for (auto iter : val)
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