#ifndef SORT_FIND_HPP
#define SORT_FIND_HPP
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <limits.h>

using namespace std;

template<typename T>
// 数组索引
int64_t findPosArray(T ar[], int n, T element)  //查找元素并返回位置下标，find(数组，长度，元素)
{
	int64_t i = 0;
	int64_t index=-1;//原始下标，没找到元素返回-1
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
int64_t findPosVector(vector<int> input, int number);

// 逆序
bool Reverse(int a, int b);

template<typename T>
// 二分查找
int64_t search_Binary_left(const vector<T>& v, T value, int64_t low = 0, int64_t high = INT64_MAX) // search_Binary_left(向量, 要找的数，左索引)
{
	int64_t vectorSize = v.size()-1;
	high = min(vectorSize, INT64_MAX);
	int64_t mid = (low + high) / 2;
	while (low <= high) {
		if (v[mid] == value) {
			return mid;
		} else if (value < v[mid]) {
			high = mid-1;
		} else {
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}

	if (high < 0) {
		high = 0;
	}
	
	return high;
}

template<typename T>
int64_t search_Binary_right(const vector<T>& v, T value,  int64_t low = 0, int64_t high = INT64_MAX)  // search_Binary_right(vector, number to find, left index)
{
	int64_t vectorSize = v.size()-1;
	high = min(vectorSize, INT64_MAX);
	int64_t mid = (low + high) / 2;
	while (low <= high) {
		if (v[mid] == value) {
			return mid;
		} else if (value < v[mid]) {
			high = mid-1;
		} else {
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}

	if (low > static_cast<int64_t>(v.size() - 1)) {
		low = (v.size() - 1);
	}
	
	return low;
}

#endif