#ifndef STRIP_SPLIT_JOIN_HPP
#define STRIP_SPLIT_JOIN_HPP
#include <string>
#include <vector>

using namespace std;


/*
	ȥ���ַ������˵������ַ�
	str -> �ַ���
	ch -> �����ַ�
*/
string strip(const string & str, char ch=' ');


/*
	���ַ������в��
	str -> �ַ���
	delim -> ��ֵ��ַ�
*/
vector<string> split(const string & str, const string & delim);


/*
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
*/
string join(vector<int> & val, string delim);


/*
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
*/
string join(vector<long int> & val, string delim);


/*
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
*/
string join(vector<float> & val, string delim);


/*
	�ַ���join�ϲ�
	val -> ��Ҫjoin������
	delim -> ����Ԫ��֮�������ַ�
*/
string join(vector<string> & val, string delim);

#endif