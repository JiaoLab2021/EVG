# include "../include/sort_find.hpp"

using namespace std;


// 向量索引
int findPosVector(vector <int> input, int number)
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
bool Reverse(int a, int b)
{
    return a > b; //升序排列，如果改为return a<b，则为降序
}
