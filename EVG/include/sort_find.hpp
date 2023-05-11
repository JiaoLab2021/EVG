#ifndef sort_find_hpp
#define sort_find_hpp
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <limits.h>

using namespace std;

template<typename T>

// 数组索引
int findPosArray(T ar[], int n, T element)//查找元素并返回位置下标，find(数组，长度，元素)
{
	int i = 0;
	int index=-1;//原始下标，没找到元素返回-1
	for (i = 0; i <n; i++)
	{
		if (element ==ar[i])
		{
			index=i;//记录元素下标
		}
	}
	return index;//返回下标
}

// 向量索引
int findPosVector(vector <int> input , int number)
{
    vector<int>::iterator iter=std::find(input.begin(),input.end(),number);//返回的是一个迭代器指针
    if(iter == input.end())
    {
        return -1;
    } else
	{
        return std::distance(input.begin(),iter);
    }
}

// 逆序
bool Reverse(int a,int b)
{
    return a > b; //升序排列，如果改为return a>b，则为降序
}

template<typename T>
// 二分查找
int search_Binary_left(vector<T>v, T value, int low = 0, int high = INT_MAX) // search_Binary_left(向量, 要找的数，左索引)
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
int search_Binary_right(vector<T>v, T value,  int low = 0, int high = INT_MAX) // search_Binary_right(向量, 要找的数，左索引)
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