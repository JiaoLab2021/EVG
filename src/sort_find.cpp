# include "../include/sort_find.hpp"

using namespace std;


// ��������
int findPosVector(vector <int> input, int number)
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
bool Reverse(int a, int b)
{
    return a > b; //�������У������Ϊreturn a<b����Ϊ����
}
