// 按照 exp(-(x-1)^2) 在 0-2之间抽样
//metropolis 

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<time.h>
#include<fstream>

using namespace std;

double w(double x) // exp(-(x-1)^2)
{
    return exp(-(x-1)*(x-1));
}

double g_01()
{
    return 1.0*rand()/RAND_MAX; 
}

double g_02() // random number gerenator
{
    return 2.0*rand()/RAND_MAX;
}

int main()
{
    srand((unsigned)time(NULL));
    const double dw=1.8;//windows 
     double N=100000;//total number of MC
     double L=10000;// previous steps before balance
    const double x0=0.99;//begin point
    int accept=0;
    double newx,xx;
    xx=x0;newx=xx;
    ofstream file1;
    file1.open("MC.dat",ios::out);
    for(int i=0;i<N;i++)
    {
        double dx=dw*(g_02()-1);
        newx=xx+dx;
        if(newx>2)  newx=newx-2;
        if(newx<0)  newx=2+newx;
        if(w(newx)/w(xx)>g_01())
        {
            xx=newx;
            accept++;
        }
        if(i>L-1)
            {
                file1<<newx<<endl;
            }
    } 
cout<<double(accept)/N;

    
 /*  for(int j=0;j<100;j++)
   {
       cout<<g_02()-1<<endl;
   }*/
    system("pause");
    return 0;
}