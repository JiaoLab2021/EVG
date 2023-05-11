#ifndef sort_find_hpp
#define sort_find_hpp
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <limits.h>

using namespace std;

template<typename T>

// ��������
int findPosArray(T ar[], int n, T element)//����Ԫ�ز�����λ���±꣬find(���飬���ȣ�Ԫ��)
{
	int i = 0;
	int index=-1;//ԭʼ�±꣬û�ҵ�Ԫ�ط���-1
	for (i = 0; i <n; i++)
	{
		if (element ==ar[i])
		{
			index=i;//��¼Ԫ���±�
		}
	}
	return index;//�����±�
}

// ��������
int findPosVector(vector <int> input , int number)
{
    vector<int>::iterator iter=std::find(input.begin(),input.end(),number);//���ص���һ��������ָ��
    if(iter == input.end())
    {
        return -1;
    } else
	{
        return std::distance(input.begin(),iter);
    }
}

// ����
bool Reverse(int a,int b)
{
    return a > b; //�������У������Ϊreturn a>b����Ϊ����
}

template<typename T>
// ���ֲ���
int search_Binary_left(vector<T>v, T value, int low = 0, int high = INT_MAX) // search_Binary_left(����, Ҫ�ҵ�����������)
{
	int vectorSize = v.size()-1;
	high = min(vectorSize, INT_MAX);
	int mid = (low + high) / 2;
	while (low <= high)
	{
		if (v[mid] == value)
		{
			return mid;
		}
		else if (value < v[mid])
		{
			high = mid-1;
		}
		else
		{
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}
	if (high < 0)
	{
		high = 0;
	}
	
	return high;
}

template<typename T>
int search_Binary_right(vector<T>v, T value,  int low = 0, int high = INT_MAX) // search_Binary_right(����, Ҫ�ҵ�����������)
{
	int vectorSize = v.size()-1;
	high = min(vectorSize, INT_MAX);
	int mid = (low + high) / 2;
	while (low <= high)
	{
		if (v[mid] == value)
		{
			return mid;
		}
		else if (value < v[mid])
		{
			high = mid-1;
		}
		else
		{
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}
	if (low > (v.size() - 1))
	{
		low = (v.size() - 1);
	}
	
	return low;
}

#endif