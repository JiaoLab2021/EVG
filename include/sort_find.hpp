#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

// Binary search in a sorted vector: index of an element equal to value, else of the last element below it (at least 0).
template<typename T>
int64_t search_Binary_left(const std::vector<T>& v, T value, int64_t low = 0, int64_t high = INT64_MAX)
{
    int64_t vectorSize = v.size()-1;
    high = std::min(vectorSize, INT64_MAX);
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

// Binary search in a sorted vector: index of an element equal to value, else of the first element above it (at most size-1).
template<typename T>
int64_t search_Binary_right(const std::vector<T>& v, T value,  int64_t low = 0, int64_t high = INT64_MAX)
{
    int64_t vectorSize = v.size()-1;
    high = std::min(vectorSize, INT64_MAX);
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
