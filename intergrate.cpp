#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;
double (*func)(double x);//函数指针 全局变量

double laginte3(double (*f)(double x),double x0,double x1) //三点拉格朗日插值 积分
{
    const int N=500;//标准去间分隔
    int n;//实际区间分割
    double S=0;
    double xx;//小区间长度
    if((x1-x0)>5)  //根据积分长度分割区间
    {
        n=(x1-x0)/5*N;
    }
    xx=(x1-x0)/n;
    for(int i=0;i<n;) //积分
    {
        S+=xx*(f(x0+xx*i)+4*f(x0+xx*(i+1))+f(x0+xx*(i+2)));
        i+=2;
    }
    return S/3.0;
} 

double function(double x)  //可以为任何其他函数
{
    double ss;
    ss=exp(x);
    return ss;
}



int main()
{
    double s,x0,x1;
    cout<<"x0 x1"<<endl;
    cin>>x0>>x1;
    s=laginte3(function,x0,x1);
    cout<<s;


    system("pause");
    return 0;
}