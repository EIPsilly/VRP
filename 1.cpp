#include<bits/stdc++.h>
using namespace std;

int main()
{
    vector<vector<int> > layer;
    int n = 10;
    layer.resize(n);
    for (int i = 0;i<n;i++) layer[i].resize(n);
    for (int i = 0;i<n;i++)
        for (int j = 0;j<n;j++)
        {
            layer[i][j]=rand();
        }
    
    layer.clear();
    cout<<layer.size()<<'\n';
    for (int i = 0;i<n;i++)
        for (int j = 0;j<n;j++)
        {
            cout<<layer[i][j];
        }
}