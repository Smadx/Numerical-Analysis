#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
namespace POL
{
    template<typename T>
    class POL
    {
        private:
            int rank;
            std::vector<T> A;
        public:
            POL();
            POL(int Rank,T a[]);
            void show();//输出表达式
            bool IsHiger(const POL &p);//比较两个多项式阶数的高低
            POL operator+(POL p);//多项式加法
            POL operator-(POL p);//多项式减法
            POL operator-();//取反
            POL operator*(T k);//多项式数乘
            POL operator*(POL p);//多项式乘法
            POL multiply(int times,T s);//多项式乘单项式
            T operator()(T x);//求值
            POL derivative(int n=1);//求导
            POL integral();//原函数
            T integral(T a,T b);//积分
    };
    template<typename T>
    POL<T>::POL()
    {
        rank=0;
        A=std::vector<T>(1);
    }
    template<typename T>
    POL<T>::POL(int Rank,T a[])
    {
        rank=Rank;
        A.resize(rank+1);
        for(int i=0;i<=rank;i++)
        {
            A[i]=a[i];
        }
    }
    template<typename T>
    void POL<T>::show()
    {
        cout<<"p(x)=";
        bool isFirst=true;
        for(int i=rank;i>=0;i--)
        {
            if(A[i]!=0)
            {
                if(isFirst) isFirst=false;
                else if (A[i]>0) cout<<'+';
                if(A[i]!=1||(i==0))
                {
                    if(A[i]==-1)
                        cout<<'-';
                    else
                        cout<<A[i];
                }    
                if (i!=0)
                {
                    cout<<'x';
                    if(i>1)
                        cout<<'^' <<i;
                }
            }
        }
        cout<<endl;
    }
    template<typename T>
    bool POL<T>::IsHiger(const POL &p)
    {
        if(p.rank>rank)
        {
            return false;
        }
        return true;
    }
    template<typename T>
    POL<T> POL<T>::operator+(POL p)
    {
        POL temp;
        temp.rank=rank;
        temp.A=A;
        POL result;
        int i;
        if(IsHiger(p))
        {
            i=rank;
            p.A.resize(i+1);
            for(int j=p.rank+1;j<=temp.rank;j++)
            {
                p.A[j]=0;
            }
        }
        else
        {
            i=p.rank;
            temp.A.resize(i+1);
            for(int j=temp.rank+1;j<=p.rank;j++)
            {
                temp.A[j]=0;
            }
        }
        result.rank=i;
        result.A.resize(i+1);
        for(int j=0;j<=i;j++)
        {
            result.A[j]=temp.A[j]+p.A[j];
        }
        int k=result.rank;
        while(result.A[k]==0)
        {
            k--;
            result.A.pop_back();
        }
        result.rank=k;
        return result;
    }
    template<typename T>
    POL<T> POL<T>::operator-(POL p)
    {
        POL result;
        int i;
        bool changed;
        if(IsHiger(p))
        {
            i=rank;
            p.A.resize(i+1);
            for(int j=p.rank+1;j<=rank;j++)
            {
                p.A[j]=0;
            }
        }
        else
        {
            i=p.rank;
            A.resize(i+1);
            for(int j=rank+1;j<=p.rank;j++)
            {
                A[j]=0;
            }
            changed=true;
        }
        result.rank=i;
        result.A.resize(i+1);
        for(int j=0;j<=i;j++)
        {
            result.A[j]=A[j]-p.A[j];
        }
        int k=result.rank+1;
        while(result.A[k]==0)
        {
            k--;
            result.A.pop_back();
        }
        result.rank=k-1;
        if(changed)
        {
            A.resize(rank+1);
        }
        return result;
    }
    template<typename T>
    POL<T> POL<T>::operator-()
    {
        for(int i=0;i<=rank;i++)
        {
            A[i]=-A[i];
        }
        return *this;
    }
    template<typename T>
    POL<T> POL<T>::operator*(T k)
    {
        POL result;
        result.rank=rank;
        result.A.resize(rank+1);
        for(int i=0;i<=rank;i++)
        {
            result.A[i]=k*A[i];
        }
        return result;
    }
    template<typename T>
    POL<T> POL<T>::multiply(int times,T s)
    {
        POL result;
        bool mul[rank+times+1]={false};
        result.rank=rank+times;
        result.A.resize(rank+times+1);
        for(int i=0;i<=rank;i++)
        {
            int t=i+times;
            mul[t]=true;
            result.A[t]=s*A[i];
        }
        for(int j=0;j<rank+times+1;j++)
        {
            if(!mul[j])
            {
                result.A[j]=0;
            }
        }
        return result;
    }
    template<typename T>
    POL<T> POL<T>::operator*(POL p)
    {
        POL result,mid;
        for(int i=0;i<=rank;i++)
        {
            mid=p.multiply(i,A[i]);
            result=result+mid;
        }
        return result;
    }
    template<typename T>
    T POL<T>::operator()(T x)
    {
        double result;
        for(int i=rank;i>=0;i--)
        {
            if(i==rank)
            {
                result=A[rank];
                continue;
            }
            result=result*x+A[i];
        }
        return result;
    }
    template<typename T>
    POL<T> POL<T>::derivative(int n)
    {
        POL result;
        result.rank=rank-n;
        result.A=A;
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<=rank-i;j++)
            {
                result.A[j]=result.A[j+1]*(j+1);
            }
            result.A.resize(rank-i);
        }
        return result;
    }
    template<typename T>
    POL<T> POL<T>::integral()
    {
        POL result;
        result.rank=rank+1;
        result.A=A;
        result.A.resize(rank+2);
        for(int i=result.rank;i>0;i--)
        {
            result.A[i]=result.A[i-1]/i;
        }
        result.A[0]=0;
        return result;
    }
    template<typename T>
    T POL<T>::integral(T a,T b)
    {
        T result;
        POL temp=integral();
        result=temp(b)-temp(a);
        return result;
    }
};