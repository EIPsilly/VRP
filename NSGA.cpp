#include<bits/stdc++.h>
using namespace std;

const int popsize = 100;// 种群大小
const int generation = 500;// 迭代次数
const int dimension = 10;// 基因个数
const int NumFun = 2;// 目标函数的个数
const int inf = 1e9;// 无穷值

// 单个个体
class individual{
    public:
        individual();
        void init();// 个体初始化
        void clear();
        void calc();// 计算目标函数的值
        friend class population;
        friend bool dominate(const individual &a,const individual &b);// 判断a是否支配b
    private:
        double gene[dimension];// 个体所携带的基因
        double fvalue[NumFun];// 函数值
        int ni;// 支配该个体的个数
        vector<int> si;// 当前个体所支配的解个体的集合
        int rank;// 支配排序分层
        double crowding_distance;// 拥挤度
};

individual::individual()
{
    memset(gene,0,sizeof(gene));
    memset(fvalue,0,sizeof(fvalue));
    ni = rank = crowding_distance = 0;
}

// 判断a是否支配b
bool dominate(const individual &a,const individual &b)
{
    bool lt = 0;
    for (int i = 0;i<NumFun;i++)
    {
        if (a.fvalue[i] > b.fvalue[i]) return false;
        if (a.fvalue[i] < b.fvalue[i]) lt = 1;
    }
    return lt;
}

void individual::init()
{
    gene[0] = (double)rand()/(double)(RAND_MAX+1.0);
    for (int i = 1;i<dimension;i++)
        gene[i] = (double)rand()/(double)(RAND_MAX) * 10.0 - 5;
}

// 计算目标函数的值
void individual::calc()
{
    int n = dimension;
    fvalue[0] = gene[0];
    double g = 91;
    for (int i = 1;i<n;i++)
        g += gene[i] * gene[i] - 10 * cos(4 * M_PI * gene[i]);
    fvalue[1] = g*(1 - sqrt(fvalue[0] / g));
}

void individual::clear()
{
    ni = rank = crowding_distance = 0;
    si.clear();
}

class population{
    public:
        void fast_nondominated_sort();//快速非支配排序
        void Elite_Strategy();// 精英策略
        void merge();// 合并P和Q成为R
        void Next_generation();// 产生下一代
        void calc_crowding_distance(int);// 计算拥挤度
        void start();// 开始迭代
        void write();
    private:
        vector<individual> P;// 父代种群
        vector<individual> Q;// 子代种群
        vector<individual> R;// 子代和父代种群
        vector<int> layer[popsize<<1];// 支配排序分层
        int rank;// 层数
};

void population::fast_nondominated_sort()
{
    // 计算每个个体与其他个体的支配关系
    for (int i = 0;i<R.size();i++)
    {
        for (int j = i+1;j<R.size();j++)
        {
            if (dominate(R[i],R[j]))
            {
                R[j].ni++;
                R[i].si.push_back(j);
            }
            else if (dominate(R[j],R[i]))
            {
                R[i].ni++;
                R[j].si.push_back(i);
            }
        }
        if (R[i].ni == 0)
            layer[0].push_back(i);
    }
    // 分级
    rank = 0;
    while (layer[rank].size() != 0)
    {
        for (vector<int>::iterator i = layer[rank].begin();i != layer[rank].end();i++)
        {
            R[*i].rank = rank;
            for (vector<int>::iterator j = R[*i].si.begin(); j != R[*i].si.end();j++)
            {
                R[*j].ni--;
                if (R[*j].ni == 0)
                    layer[rank+1].push_back(*j);
            }
        }
        rank++;
    }
}

// 计算每一层的拥挤度
void population::calc_crowding_distance(int index)
{
    for (int j = 0;j<NumFun;j++)
    {
        sort(layer[index].begin(),layer[index].end(),[&](const int a,const int b){
            return R[a].fvalue[j] < R[b].fvalue[j];
        });
        // 边界的两个个体拥挤度为无穷

        for (vector<int>::iterator it = layer[index].begin();it != layer[index].end();it++)
        {
            if (it == layer[index].begin() || (it + 1) == layer[index].end())
                R[*it].crowding_distance = inf;
            else
                R[*it].crowding_distance += fabs(R[*(it+1)].fvalue[j] - R[*(it-1)].fvalue[j]);
        }
    }    
}

// 精英策略
void population::Elite_Strategy()
{
    for (int i = 0;i<rank;i++)
    {
        calc_crowding_distance(i);
        if (P.size() + layer[i].size() <= popsize)
        {
            for (vector<int>::iterator it = layer[i].begin();it != layer[i].end();it++)
                P.push_back(R[*it]);
        }
        else
        {
            // layer中的标号按照拥挤度从大到小排序
            sort(layer[i].begin(),layer[i].end(),[&](const int a,const int b){return R[a].crowding_distance > R[b].crowding_distance;});

            for (vector<int>::iterator it = layer[i].begin();it != layer[i].end() && P.size() < popsize;it++)
                P.push_back(R[*it]);
            break;
        }
    }
    for (vector<individual>::iterator it = P.begin();it != P.end();it++)
        it->clear();
}

void population::merge()
{
    R.clear();
    for (vector<individual>::iterator it = P.begin();it != P.end();it++)
    {
        it->clear();
        R.push_back(*it);
    }
    P.clear();

    for (vector<individual>::iterator it = Q.begin();it != Q.end();it++)
    {
        it->clear();
        R.push_back(*it);
    }
    Q.clear();

    for (int i = 0;i<rank;i++) layer[i].clear();
    rank = 0;
}

void population::Next_generation()
{
    vector<individual> tmp;
    // 二元竞赛选择
    for (int i = 0;i<popsize;i++)
    {
        int x = rand() % popsize,y = rand() % popsize;
        while (y == x) y = rand() % popsize;

        if (P[x].rank < P[y].rank) tmp.push_back(P[x]);
        else if (P[y].rank < P[y].rank) tmp.push_back(P[y]);
        else if (P[x].crowding_distance > P[y].crowding_distance) tmp.push_back(P[x]);
        else tmp.push_back(P[y]);
    }

    P.swap(tmp);

    // 模拟二进制交叉
    for (int i = 0; i < popsize / 2; i++) 
    {
        double probability = (double)rand() / (double)(RAND_MAX);
        if (probability <= 0.8)
        {
            double u = (double)rand() / (double)(RAND_MAX + 1),b;
            if (u <= 0.5) b = pow(2*u,1.0/21.0);
            else b = 1.0 / pow(2.0 - 2 * u,1.0/21.0);

            int x = rand() % popsize,y = rand() % popsize;
            while (y == x) y = rand() % popsize;

            individual c1,c2;
            for (int j = 0;j<dimension;j++)
            {
                c1.gene[j] = (P[x].gene[j] + P[y].gene[j]) / 2.0 - b * (P[y].gene[j] - P[x].gene[j]) / 2.0;
                c2.gene[j] = (P[x].gene[j] + P[y].gene[j]) / 2.0 + b * (P[y].gene[j] - P[x].gene[j]) / 2.0;
                c1.gene[j] = max(0.0,c1.gene[j]);
                c1.gene[j] = min(1.0,c1.gene[j]);
                c2.gene[j] = max(0.0,c2.gene[j]);
                c2.gene[j] = min(1.0,c2.gene[j]);
            }
            Q.push_back(c1);
            Q.push_back(c2);
        }
    }

    // 多项式变异
    for (int i = 0;i<popsize;i++)
    {
        double probability = (double)rand() / (double)(RAND_MAX);
        if (probability >= 0.95)
        {
            individual c1;
            for(int j=0;j<dimension;j++)
            {
                int x = rand() % popsize;
                double u = (double)rand() / (double)(RAND_MAX + 1),b;
                if(u<0.5)
                    b = pow(2*u,1.0/21.0) - 1;
                else
                    b = 1 - pow(2-2*u,1.0/21.0);

                c1.gene[j]=P[x].gene[j]+(1.0-0.0)*u;
                c1.gene[j] = max(0.0,c1.gene[j]);
                c1.gene[j] = min(1.0,c1.gene[j]);
            }
            Q.push_back(c1);
        }
    }
    // 计算目标函数值
    for (vector<individual>::iterator it = Q.begin();it != Q.end();it++)
        it->calc();
}

void population::start()
{
    for (int i = 0;i<popsize;i++)
    {
        individual tmp;
        tmp.init();
        tmp.calc();
        R.push_back(tmp);
    }
    rank = 0;
    for (int times = 0;times < generation;times++)
    {
        // cout<<"generations "<<times<<'\n';
        // 非支配排序
        fast_nondominated_sort();
        // 精英策略从R中筛选popsize个放入P中
        Elite_Strategy();
        // 生成子代Q
        Next_generation();
        // 合并P和Q成为R
        merge();
    }
}

void population::write()
{
    freopen("./ZDT4/my.txt","w",stdout);
    for (int i = 0;i<popsize;i++)
        cout<<P[i].fvalue[0]<<' '<<P[i].fvalue[1]<<'\n';
    fclose(stdout);
}

int main(){
    clock_t st,ed;
    srand(time(NULL));
    st = clock();
    population solve;
    solve.start();
    ed = clock();
    cout<<(double)(ed - st) / CLOCKS_PER_SEC;
    solve.write();
    return 0;
}