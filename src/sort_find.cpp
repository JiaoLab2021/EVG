# include "../include/sort_find.hpp"

using namespace std;


// Vector index
int64_t findPosVector(vector <int> input, int number)
{
    vector<int>::iterator iter=std::find(input.begin(),input.end(),number);  // Returns an iterator pointer
    if(iter == input.end())
    {
        return -1;
    } else
	{
        return std::distance(input.begin(),iter);
    }
}

// Reverse sequence
bool Reverse(int a, int b)
{
    return a > b; // In ascending order, if return a<b is changed, it is in descending order
}
