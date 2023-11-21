#include "Matrix.h"
using namespace std;
int main()
{
    //创建一个类
    double a[3][3]={{1,2,3},{4,5,6},{7,8,9}};
    Matrix::Matrix<double> A(3,3,(double **)a);
    cout<<"A:"<<A<<endl;
    double b[3][3]={{-1,2,3},{4,-5,6},{7,8,-9}};
    Matrix::Matrix<double> B(3,3,(double **)b);
    cout<<"B:"<<B<<endl;
    double c[3]={1,2,3};
    Vector::Vector<double> C(3,c);
    cout<<"C:"<<C<<endl;
    //使用类的成员函数
    /*访问类的成员*/
    //访问矩阵的元素
    cout<<A[1][1]<<endl;
    //矩阵的行数和列数
    cout<<A.row()<<endl;
    cout<<A.col()<<endl;
    //矩阵的转置
    cout<<A.Trans()<<endl;
    //矩阵的逆
    cout<<B.Inverse()<<endl;
    //矩阵的行列式
    cout<<B.Det()<<endl;
    /*使用类的重载运算符*/
    //矩阵的加减法
    cout<<A+B<<endl;
    cout<<A-B<<endl;
    //矩阵的数乘
    cout<<A*2<<endl;
    //矩阵的乘法
    cout<<A*B<<endl;
    //向量与矩阵相乘
    cout<<C*A<<endl;
    return 0;
}