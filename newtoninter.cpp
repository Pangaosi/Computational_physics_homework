#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

double newton(double x[],double y[],int n,double xx) //n个点,xx
{   double yx=y[0]; //判断初始
    if(n==0)
    {
        return yx;
    }
    else 
    {
        double  oi,Nii;// oi 定义见下；Nii为 N_(n-1);
       Nii=newton(x,y,n-1,x[n]);
       oi=y[n]-Nii;
       for(int j=0;j<n;j++)
       {
           oi*=(xx-x[j])/(x[n]-x[j]);
       }
       yx=newton(x,y,n-1,xx)+oi;  //迭代

        return yx;
    }

}

int main()
{
    const int N=10;  // 10个点，原始数据
    double tt=4*atan(1)/N;
    double x[N+1],y[N+1];
    for(int i=0;i<N+1;i++)
    {
        x[i]=i*tt;
        y[i]=cos(x[i]);
    }
    int n=N*100;  //插值点
    double ttt=tt/100;
    double yx[n];
    ofstream file1;
    file1.open("newtoninter.dat",ios::out);
    for(int j=0;j<n;j++)
    {
        yx[j]=newton(x,y,N,(j+1)*ttt);
        file1<<(j+1)*ttt<<'\t'<<yx[j]<<endl;
    }
    cout<<"success"<<endl;
    file1.close();
    cout<<newton(x,y,N,atan(1));
    system("pause");
    return 0;

}