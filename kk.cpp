#include<iostream>
#include<fstream>
#include<cmath>
 #define _USE_MATH_DEFINES

using namespace std;

double interpolation(double x[],double y[],int n,double xx)
{   
     double AA=1;
     double S=0;
     for(int i=0;i<n;i++)
     {
         for(int j=0;j<n;j++)
         {
            if(j!=i)
            {
                AA*=(xx-x[j])/(x[i]-x[j]);
            }
         }
         S+=AA*y[i];AA=1;
     }
     return S;
}

int main()
{
    const double t=atan(1)/2.0;
    double x[5],y[5];
    double m[1000],n[1000];
    double v[1000],w[1000];
    for(int k=0;k<5;k++)
    {
        x[k]=k*t;
        y[k]=cos(x[k]);
    }
    for(int p=0;p<=1000;p++)
    {  double tt=2*atan(1)/1001.0;
        m[p]=p*tt;n[p]=cos(m[p]);
        v[p]=m[p];w[p]=interpolation(x,y,5,v[p]);
    }
    ofstream file1;ofstream file2;
    file1.open("d1.dat",ios::out);
    file2.open("d2.dat",ios::out);
    double ss=0;
    for(int q=0;q<1000;q++)
    {
        file1<<m[q]<<'\t'<<n[q]<<endl;
        file2<<v[q]<<'\t'<<w[q]<<endl;
        ss+=0.001;
    }
    file1.close();
    file2.close();
    cout<<ss;
    
    
    system("pause");
    return 0;
}
