#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;
namespace COMPLEX
{
    class Complex
    {
        private:
            double real;
            double imag;
        public:
            Complex();//构造0
            Complex(double x,double y,bool angel=0);//构造x+yi,angel=1时y为辐角
            double Re();//获取实部
            double Im();//获取虚部
            double mod();//获取模
            double arg();//获取辐角
            Complex operator+(Complex B);//加法
            Complex operator-(Complex B);//减法
            Complex operator-();//取反
            Complex operator*(Complex B);//乘法
            Complex operator/(Complex B);//除法
            template<typename t>
            friend Complex operator*(Complex A,t lambda);//实数数乘(友元)
            template<typename t>
            friend Complex operator*(t lambda,Complex A);//实数数乘(友元)
            template<typename t>
            friend Complex operator/(Complex A,t lambda);//实数数除(友元)
            template<typename t>
            friend Complex operator/(t lambda,Complex A);//实数数除(友元)
            template<typename t>
            friend Complex operator^(Complex A,t lambda);//幂(友元)
            template<typename t>
            friend Complex operator^(t lambda,Complex A);//幂(友元)
            friend std::ostream& operator<<(std::ostream& os,const Complex &A);//输出(友元)
            friend std::ofstream& operator<<(std::ofstream& fout,const Complex &A);//输出(友元)
    };
    Complex::Complex()
    {
        real=0;
        imag=0;
    }
    Complex::Complex(double x,double y,bool angel)
    {
        if(angel==0)
        {
            real=x;
            imag=y;
        }
        else
        {
            real=x*cos(y);
            imag=x*sin(y);
        }
    }
    double Complex::Re()
    {
        return real;
    }
    double Complex::Im()
    {
        return imag;
    }
    double Complex::mod()
    {
        return sqrt(real*real+imag*imag);
    }
    double Complex::arg()
    {
        return atan2(imag,real);
    }
    Complex Complex::operator+(Complex B)
    {
        return Complex(real+B.real,imag+B.imag);
    }
    Complex Complex::operator-(Complex B)
    {
        return Complex(real-B.real,imag-B.imag);
    }
    Complex Complex::operator-()
    {
        return Complex(-real,-imag);
    }
    Complex Complex::operator*(Complex B)
    {
        return Complex(real*B.real-imag*B.imag,real*B.imag+imag*B.real);
    }
    Complex Complex::operator/(Complex B)
    {
        return Complex((real*B.real+imag*B.imag)/(B.real*B.real+B.imag*B.imag),(imag*B.real-real*B.imag)/(B.real*B.real+B.imag*B.imag));
    }
    template<typename t>
    Complex operator*(Complex A,t lambda)
    {
        return Complex(A.real*lambda,A.imag*lambda);
    }
    template<typename t>
    Complex operator*(t lambda,Complex A)
    {
        return Complex(A.real*lambda,A.imag*lambda);
    }
    template<typename t>
    Complex operator/(Complex A,t lambda)
    {
        return Complex(A.real/lambda,A.imag/lambda);
    }
    template<typename t>
    Complex operator/(t lambda,Complex A)
    {
        return Complex(A.real/lambda,A.imag/lambda);
    }
    template<typename t>
    Complex operator^(Complex A,t lambda)
    {
        return Complex(pow(A.mod(),lambda)*cos(lambda*A.arg()),pow(A.mod(),lambda)*sin(lambda*A.arg()));
    }
    template<typename t>
    Complex operator^(t lambda,Complex A)
    {
        return Complex(pow(A.mod(),lambda)*cos(lambda*A.arg()),pow(A.mod(),lambda)*sin(lambda*A.arg()));
    }
    std::ostream& operator<<(std::ostream& os,const Complex &A)
    {
        if(A.imag*A.imag+A.real*A.real<1e-10)
        {
            os<<0;
            return os;
        }
        if(A.imag>0)
            os<<A.real<<"+"<<A.imag<<"i";
        else if(A.imag==0)
            os<<A.real;
        else
            os<<A.real<<A.imag<<"i";
        return os;
    }
    std::ofstream& operator<<(std::ofstream& fout,const Complex &A)
    {
        if(A.imag*A.imag+A.real*A.real<1e-10)
        {
            fout<<0;
            return fout;
        }
        if(A.imag>0)
            fout<<A.real<<"+"<<A.imag<<"i";
        else if(A.imag==0)
            fout<<A.real;
        else
            fout<<A.real<<A.imag<<"i";
        return fout;
    }
}