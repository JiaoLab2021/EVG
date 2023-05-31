#include "../include/strip_split_join.hpp"

using namespace std;


/*
	ȥ���ַ������˵������ַ�
	str -> �ַ���
	ch -> �����ַ�
*/
string strip(const string & str, char ch)
{
	int i = 0;
	while (str[i] == ch)  // ͷ��ch�ַ������� i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch)  // ��βch�ַ������� str.size() - 1 - j
		j--;		
	return str.substr(i, j+1-i);
}


/*
	���ַ������в��
	str -> �ַ���
	delim -> ��ֵ��ַ�
*/
vector<string> split(const string & str, const string & delim)
{
	vector<string> res;  // ���ָ������ַ����洢��vector��
	if("" == str) return  res;  // ���Ϊ�շ��ؿ�vector
	
	string strs = str + delim;  // ��չ�ַ����Է���������һ���ָ������ַ���
	size_t pos;
	size_t size = strs.size();
 
	for (int i = 0; i < size; ++i) {
		pos = strs.find(delim, i); // posΪ�ָ�����һ�γ��ֵ�λ�ã���i��pos֮ǰ���ַ����Ƿָ��������ַ���
		if( pos < size)  // ������ҵ������û�в��ҵ��ָ�����posΪstring::npos
		{
			string s = strs.substr(i, pos - i);  // ��i��ʼ����Ϊpos-i�����ַ���
			res.push_back(s);  // ���������ո�֮���и�����ַ���Ϊ���ַ���������û���ж�s�Ƿ�Ϊ�գ��������Ľ�����п��ַ��������
			i = pos + delim.size() - 1;
		}
	}
	return res;	
}


/*
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
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
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
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
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
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
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
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