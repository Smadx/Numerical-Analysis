#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;
namespace Matrix
{
    template<typename T>
    class Matrix
    {
        private:
            int x;
            int y;
            std::vector<vector<T> > M;     
        public:
            Matrix();//默认构造
            Matrix(int xx,int yy);//构造xx行yy列矩阵(空)
            Matrix(int xx,int yy,T **p);//构造xx行yy列矩阵
            Matrix(int n);//构造n阶单位阵
            void resize(int xx,int yy);//更改大小
            void set(int xx,int yy,T p);//更改矩阵的元素
            T row();//行
            T col();//列
            void show();//显示矩阵
            void elem_exchange(int a,int b,bool r=1);//初等变换:交换矩阵的第a行和第b行(r=1表示行变换)
            void elem_scalar_mult(double lambda,int a,bool r=1);//初等变换:把第a行数乘lambda(r=1表示行变换)
            void elem_add(int a,double lambda,int b,bool r=1);//初等变换:把第a行乘上lambda后加到第b行(r=1表示行变换)
            Matrix Gauss();//Gauss消元
            T Det();//行列式
            Matrix Inverse();//求逆
            Matrix Trans();//转置
            Matrix Diag();//取对角元
            Matrix loc(int a,bool r=1);//取出第a行(r=1表示取出行)
            vector<T>& operator[](int a);//重载[]
            Matrix operator+(Matrix B);//矩阵相加
            bool operator==(Matrix B);//重载等号
            Matrix operator-(Matrix B);//矩阵相减
            Matrix operator-();//取反
            Matrix operator*(T m);//数乘
            Matrix operator*(Matrix B);//矩阵乘法
            Matrix operator>>(Matrix B);//拼块
            template<typename t>
            friend Matrix<t> operator*(double m,Matrix<t> A);//数乘(友元)
            template<typename t>
            friend Matrix<t> operator*(Matrix<t> m,Matrix<t> A);//矩阵乘法(友元)
            template<typename t>
            friend Matrix<t> operator<<(Matrix<t> A,Matrix<t> B);//拼块(友元)
            template<typename t>
            friend std::ostream &operator<<(std::ostream &os,const Matrix<t> &A);//输出到cmd
            template<typename t>
            friend std::ofstream &operator<<(std::ofstream &of,const Matrix<t> &A);//输出到文件
    };
    template<typename T>
    Matrix<T>::Matrix()
    {
        x=0;
        y=0;
        M=vector<vector<double> >(0);
    }
    template<typename T>
    Matrix<T>::Matrix(int xx,int yy)
    {
        x=xx;
        y=yy;
        M.resize(xx);
        for(int i=0;i<xx;i++)
        {
            M[i].resize(yy);
        }
    }
    template<typename T>
    Matrix<T>::Matrix(int xx,int yy,T **p)
    {
        x=xx;
        y=yy;
        for(int i=0;i<x;i++)
        {
            M.resize(x);
            for(int j=0;j<y;j++)
            {
                M[i].resize(y);
                M[i][j]=*((T *)p+i*y+j);
            }
        }
    }
    template<typename T>
    Matrix<T>::Matrix(int n)
    {
        x=n;
        y=n;
        M.resize(n);
        for(int i=0;i<n;i++)
        {
            M[i].resize(n);
            for(int j=0;j<n;j++)
            {
                if(i==j)
                    M[i][j]=1;
            }
        }
    }
    template<typename T>
    void Matrix<T>::resize(int xx,int yy)
    {
        x=xx;
        y=yy;
        M.resize(xx);
        for(int i=0;i<xx;i++)
        {
            M[i].resize(yy);
        }
    }
    template<typename T>
    void Matrix<T>::set(int xx,int yy,T p)
    {
        M[xx][yy]=p;
    }
    template<typename T>
    T Matrix<T>::row()
    {
        return x;
    }
    template<typename T>
    T Matrix<T>::col()
    {
        return y;
    }
    template<typename T>
    void Matrix<T>::show()
    {
        for(int i=0;i<x;i++)
        {
            for(int j=0;j<y;j++)
            {
                cout<<M[i][j]<<"  ";
            }
            cout<<endl;
        }
    }
    template<typename T>
    void Matrix<T>::elem_exchange(int a,int b,bool r)//交换行/列
    {
        if(r==1)//r=1表示行变换,否则为列变换
        {
            vector<T> temp=M[a-1];
            M[a-1]=M[b-1];
            M[b-1]=temp;
        }
        else
        {
            vector<T> temp;
            temp.resize(x);
            for(int i=0;i<x;i++)
            {
                temp[i]=M[i][a-1];
                M[i][a-1]=M[i][b-1];
                M[i][b-1]=temp[i];
            }
        }
    }
    template<typename T>
    void Matrix<T>::elem_scalar_mult(double lambda,int a,bool r)//某一行/列数乘
    {
        if(r==1)
        {
            for(int i=0;i<y;i++)
            {
                M[a-1][i]=lambda*M[a-1][i];
            }
        }
        else
        {
            for(int i=0;i<x;i++)
            {
                M[i][a-1]=lambda*M[i][a-1];
            }
        }
    }
    template<typename T>
    void Matrix<T>::elem_add(int a,double lambda,int b,bool r)
    {
        if(r==1)
        {
            for(int i=0;i<y;i++)
            {
                M[b-1][i]+=lambda*M[a-1][i];
            }
        }
        else
        {
            for(int i=0;i<x;i++)
            {
                M[i][b-1]+=lambda*M[i][a-1];
            }
        }
    }
    template<typename T>
    Matrix<T> Matrix<T>::Trans()
    {
        Matrix result;
        result.x=y;
        result.y=x;
        result.M.resize(y);
        for(int i=0;i<y;i++)
        {
            result.M[i].resize(x);
            for(int j=0;j<x;j++)
            {
                result.M[i][j]=M[j][i];
            }
        }
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::Diag()
    {
        Matrix result(x,y);
        result.x=x;
        result.y=y;
        for(int i=0;i<x;i++)
        {
            for(int j=0;j<y;j++)
            {
                if(i==j)
                    result.M[i][j]=M[i][j];
                else
                    result.M[i][j]=0;
            }
        }
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::loc(int a,bool r)
    {
        if(r==1)
        {
            Matrix result(1,y);
            result.M[0]=M[a];
            return result;
        }
        else
        {
            Matrix result(x,1);
            for(int i=0;i<this->row();i++)
            {
                result.M[i][0]=M[i][a];
            }
            return result;
        }
    }
    template<typename T>
    vector<T>& Matrix<T>::operator[](int a)
    {
        return M[a];
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator+(Matrix B)
    {
        Matrix result;
        if((B.x!=x)||(B.y!=y))
        {
            cout<<"add unmatched"<<endl;
            exit(1);
        }
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(y);
            for(int j=0;j<y;j++)
            {
                result.M[i][j]=B.M[i][j]+M[i][j];
            }
        }
        result.x=x;
        result.y=y;
        return result;
    }
    template<typename T>
    bool Matrix<T>::operator==(Matrix B)
    {
        const double eps = 1e-6;
        Matrix result;
        if((B.x!=x)||(B.y!=y))
        {
            cout<<"is_equal unmatched"<<endl;
            exit(1);
        }
        for(int i=0;i<x;i++)
        {
            for(int j=0;j<y;j++)
            {
                if(abs(B.M[i][j]-M[i][j])>eps) return false;
            }
        }
        return true;
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator-(Matrix B)
    {
        Matrix result;
        if((B.x!=x)||(B.y!=y))
        {
            cout<<"unmatched"<<endl;
            exit(1);
        }
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(y);
            for(int j=0;j<y;j++)
            {
                result.M[i][j]=M[i][j]-B.M[i][j];
            }
        }
        result.x=x;
        result.y=y;
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator-()
    {
        for(int i=0;i<x;i++)
        {
            for(int j=0;j<y;j++)
            {
                M[i][j]=-M[i][j];
            }
        }
        return *this;
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator*(T m)
    {
        Matrix result;
        result.x=x;
        result.y=y;
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(y);
            for(int j=0;j<y;j++)
            {
                result.M[i][j]=m*M[i][j];
            }
       }
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator*(Matrix B)
    {
        Matrix result;
        T res;
        if((y!=B.x))
        {
            cout<<"multiply unmatched"<<endl;
            exit(1);
        }
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(B.y);
            for(int j=0;j<B.y;j++)
            {
                res=0;
                for(int m=0;m<y;m++)
                {
                    res+=M[i][m]*B.M[m][j];
                }
                result.M[i][j]=res;
            }
        }
        result.x=x;
        result.y=B.y;
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::operator>>(Matrix B)
    {
        if(B.x!=x)
        {
            cout<<"unmatched"<<endl;
            exit(1);
        }
        Matrix result;
        result.x=x;
        result.y=y+B.y;
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(result.y);
            for(int j=0;j<result.y;j++)
            {
                if(j<y)
                {
                    result.M[i][j]=M[i][j];
                }
                else
                {
                    result.M[i][j]=B.M[i][j-y];
                }
            }
        }
        return result;
    }
    template<typename t>
    Matrix<t> operator*(double  m,Matrix<t> A)
    {
        return A*m;
    }
    template<typename t>
    Matrix<t> operator*(Matrix<t> m,Matrix<t> A)
    {
        return A*m;
    }
    template<typename T>
    Matrix<T> operator<<(Matrix<T> A,Matrix<T> B)
    {
        if(B.x!=A.x)
        {
            cout<<"unmatched"<<endl;
            exit(1);
        }
        Matrix<T> result;
        result.x=B.x;
        result.y=B.y+A.y;
        result.M.resize(B.x);
        for(int i=0;i<result.x;i++)
        {
            result.M[i].resize(result.y);
            for(int j=0;j<result.y;j++)
            {
                if(j<A.y)
                {
                    result.M[i][j]=A.M[i][j];
                }
                else
                {
                    result.M[i][j]=B.M[i][j-A.y];
                }
            }
        }
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::Gauss()
    {
        Matrix result=*this;
        for(int i=0;i<x-1;i++)
        {
            for(int j=i+1;j<x;j++)
            {    
                result.elem_add(i+1,-result.M[j][i]/result.M[i][i],j+1);
            }
        }
        return result;
    }
    template<typename T>
    T Matrix<T>::Det()
    {
        Matrix A=Gauss();
        T result=A.M[0][0];
        for(int i=1;i<x;i++)
        {
            result*=A.M[i][i];
        }
        return result;
    }
    template<typename T>
    Matrix<T> Matrix<T>::Inverse()
    {
        if(!(this->Det()))
        {
            cout<<"It cant be inversed"<<endl;
        }
        Matrix I(x);
        Matrix result;
        result.x=x;
        result.y=y;
        Matrix temp=*this>>I;
        temp=temp.Gauss();
        for(int i=x;i>1;i--)
        {
            for(int j=i-1;j>0;j--)
            {    
                temp.elem_add(i,-temp.M[j-1][i-1]/temp.M[i-1][i-1],j);
            }
        }
        for(int i=0;i<x;i++)
        {
            if(temp.M[i][i]!=1)
                temp.elem_scalar_mult(1/temp.M[i][i],i+1);
        }
        result.M.resize(x);
        for(int i=0;i<x;i++)
        {
            result.M[i].resize(y);
            for(int j=0;j<y;j++)
            {
                result.M[i][j]=temp.M[i][j+x];
            }
        }
        return result;
    }
    template<typename T>
    std::ostream &operator<<(std::ostream &os,const Matrix<T> &A)
    {
        for(int i=0;i<A.x;i++)
        {
            for(int j=0;j<A.y;j++)
            {
                if(fabs(A.M[i][j])<1e-6)
                    os<<0<<"  ";
                else
                    os<<A.M[i][j]<<"  ";
            }
            os<<endl;
        }
        return os;
    }
    template<typename T>
    std::ofstream &operator<<(ofstream &fout,const Matrix<T> &A)
    {
        for(int i=0;i<A.x;i++)
        {
            for(int j=0;j<A.y;j++)
            {
                if(fabs(A.M[i][j])<1e-6)
                    fout<<0<<"  ";
                else
                    fout<<A.M[i][j]<<"  ";
            }
            fout<<endl;
        }
        return fout;
    }
}
namespace Vector
{
    template<typename T>
    class Vector : public Matrix::Matrix<T>//默认为列向量
    {
    public:
        Vector():Matrix::Matrix<T>(){};//默认构造函数
        Vector(int n):Matrix::Matrix<T>(n,1){};//构造n维的零向量
        Vector(int n,T *a):Matrix::Matrix<T>(n,1,(T **)a){};//构造n维的向量，元素为a
        Vector(Matrix::Matrix<T> A);//将矩阵A转化为向量
        void resize(int n);//将向量的维数变为n
        void push_back(T a);//在向量的末尾添加元素a
        int dim();//返回向量的维数
        T norm(int p=2);//返回向量的p范数
        Vector<T> unit();//返回向量的单位向量
        T& operator[](int i);//返回向量的第i个元素
        Vector<T> operator+(Vector<T> B);//向量加法
        Vector<T> operator-(Vector<T> B);//向量减法
        Vector<T> operator-();//向量取负
        bool operator==(Vector<T> B);//向量比较
        T operator*(Vector<T> B);//向量点乘
        Vector<T> operator^(Vector<T> B);//向量叉乘
        template<typename t1,typename t2>
        friend Vector<t1> operator*(Vector<t1> A,t2 m);//向量数乘与矩阵乘法(友元)
        template<typename t1,typename t2>
        friend Vector<t1> operator*(t2 m,Vector<t1> A);//向量数乘与矩阵乘法(友元)
        template<typename t>
        friend std::ostream &operator<<(std::ostream &os,const Vector<t> &A);//向量输出
        template<typename t>
        friend std::ofstream &operator<<(std::ofstream &fout,const Vector<t> &A);//向量输出(文件)
    };
    template<typename T>
    Vector<T>::Vector(Matrix::Matrix<T> A)
    {
        if(A.row()==1)
        {
            this->resize(A.col());
            for(int i=0;i<this->row();i++)
            {
                this->set(i,0,A[0][i]);
            }
        }
        else if((A.row()!=1)&&(A.col()==1))
        {
            this->resize(A.row());
            for(int i=0;i<this->row();i++)
            {
                this->set(i,0,A[i][0]);
            }
        }
        else
        {
            cout<<"It is not a vector"<<endl;
        }
    };
    template<typename T>
    void Vector<T>::resize(int n)
    {
        Matrix::Matrix<T>::resize(n,1);
    }
    template<typename T>
    void Vector<T>::push_back(T a)
    {
        this->resize(this->dim()+1);
        this->M[this->dim()-1][0]=a;
    }
    template<typename T>
    int Vector<T>::dim()
    {
        return this->row();
    }
    template<typename T>
    T Vector<T>::norm(int p)
    {
        T result=0;
        for(int i=0;i<this->dim();i++)
        {
            result+=pow(fabs(this->M[i][0]),p);
        }
        return pow(result,1.0/p);
    }
    template<typename T>
    Vector<T> Vector<T>::unit()
    {
        return *this/this->norm();
    }
    template<typename T>
    T& Vector<T>::operator[](int i)
    {
        return this->M[i][0];
    }
    template<typename T>
    Vector<T> Vector<T>::operator+(Vector<T> B)
    {
        return *this+static_cast<Matrix::Matrix<T> >(B);
    }
    template<typename T>
    Vector<T> Vector<T>::operator-(Vector<T> B)
    {
        return *this-static_cast<Matrix::Matrix<T> >(B);
    }
    template<typename T>
    Vector<T> Vector<T>::operator-()
    {
        return -static_cast<Matrix::Matrix<T> >(*this);
    }
    template<typename T>
    bool Vector<T>::operator==(Vector<T> B)
    {
        return *this==static_cast<Matrix::Matrix<T> >(B);
    }
    template<typename T>
    T Vector<T>::operator*(Vector<T> B)
    {
        double result=0;
        for(int i=0;i<this->dim();i++)
        {
            result+=this->M[i][0]*B[i];
        }
        return result;
    }
    template<typename T>
    Vector<T> Vector<T>::operator^(Vector<T> B)
    {
        //判断维数是否为3
        if(this->dim()!=3||B.dim()!=3)
        {
            cout<<"Error: The dimension of the vector is not 3."<<endl;
            exit(0);
        }
        Vector<T> result(3);
        result[0]=this->M[1][0]*B[2]-this->M[2][0]*B[1];
        result[1]=this->M[2][0]*B[0]-this->M[0][0]*B[2];
        result[2]=this->M[0][0]*B[1]-this->M[1][0]*B[0];
        return result;
    }
    template<typename t1,typename t2>
    Vector<t1> operator*(Vector<t1> A,t2 m)
    {
        return Vector<t1>(static_cast<Matrix::Matrix<t1> >(A)*m);
    }
    template<typename t1>
    Vector<t1> operator*(Vector<t1> A,Matrix::Matrix<t1> B)
    {
        return Vector<t1>((static_cast<Matrix::Matrix<t1> >(A).Trans()) * B);
    }
    template<typename t1,typename t2>
    Vector<t1> operator*(t2 m,Vector<t1> A)
    {
        return Vector<t1>(static_cast<Matrix::Matrix<t1> >(A)*m);
    }
    template<typename t1>
    Vector<t1> operator*(Matrix::Matrix<t1> A,Vector<t1> B)
    {
        return Vector<t1>(A * static_cast<Matrix::Matrix<t1> >(B));
    }
    template<typename t>
    std::ostream &operator<<(std::ostream &os,const Vector<t> &A)
    {
        os<<static_cast<Matrix::Matrix<t> >(A).Trans();
        return os;
    }
    template<typename t>
    std::ofstream &operator<<(std::ofstream &fout,const Vector<t> &A)
    {
        fout<<static_cast<Matrix::Matrix<t> >(A).Trans();
        return fout;
    }
}