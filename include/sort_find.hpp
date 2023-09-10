#ifndef SORT_FIND_HPP
#define SORT_FIND_HPP
#include <iostream>
#include <algorithm>
#include <vector>
#include <regex>
#include <limits.h>

using namespace std;

template<typename T>
// Array index
int64_t findPosArray(T ar[], int n, T element)  // Finds the element and returns the position subscript, find(array, length, element)
{
	int64_t i = 0;
	int64_t index=-1;  // Original subscript, no element found returns -1
	for (i = 0; i <n; i++)
	{
		if (element ==ar[i])
		{
			index=i;  // Record element subscript
		}
	}
	return index;  // Returns the subscript
}

// Vector index
int64_t findPosVector(vector<int> input, int number);

// Reverse order
bool Reverse(int a, int b);

template<typename T>
int64_t search_Binary_left(const vector<T>& v, T value, int64_t low = 0, int64_t high = INT64_MAX)  // search_Binary_left(vector, number to find, left index)
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