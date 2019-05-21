#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<ctime>

using namespace std;

const int MCS=10000;//总蒙卡步数 
//const int Mpe=2000;//前期蒙卡步数
const double dw=1.8;//窗口
const double T0=150;
const double dt=0.3;

const double money[10][3]={{100,70,50},{80,65,40},{120,100,80},{90,70,50},{100,75,60},{150,120,100},{60,50,40},{85,70,55},{110,90,60},{130,100,60}};
const double mue[10]={19.920925,19.911963,9.7847135,15.002913,24.9486635,15.902663,18.081303,12.311147,15.1194365,22.040846};
//const double sigma[10]={1.814153607,1.946467664,1.935551048,1.873943848,1.849842189,2.190870125,2.096132932,1.953823092,2.050788227,1.986059336};
double X[10];//={19.1349,21.2705,9.8569,17.2835,26.5164,17.4319,17.8544,13.062,16.6191,21.5075};
double X_bestnow[10];
double profit=0;
double profit_best=0;


double costs(double x[10])
{
    double sum=0;
    for(int i=0;i<10;i++)
    {
        double dx=fabs(x[i]-mue[i])/mue[i];
        if(dx<0.05)
        sum+=money[i][0];
        else if((dx>0.05)&&(dx<0.1))
        sum+=money[i][1];
        else 
        sum+=money[i][2];
    }
    return sum;
}

double product(double quality) //返回产品售价
{
    if(quality>182&&quality<186)
    return 1500;
    else if((quality>178&&quality<182)||(quality>186&&quality<190))
    return 1200;
    else if((quality>165&&quality<178)||(quality>190&&quality<197))
    return 800;
    else 
    return 600;
}

double f(double x[10]) ////////////////////////////函数关系
{
    return  -10.8323+1.0041*x[0]+1.0043*x[1]+0.9961*x[2]+1.0013*x[3]+1.4926*x[4]+0.9958*x[5]+0.9949*x[6]+1.0010*x[7]+0.9993*x[8]+1.4403*x[9];

}

void initial(double x[10]) ////////////////////初始化
{
    for(int i=0;i<10;i++)
    {
        x[i]=mue[i]+4.0*(2.0*rand()/RAND_MAX-1);
    }
    profit=(product(f(x))-costs(x))/costs(x);
}

void change(double x[10],double y[10])
{
    for(int i=0;i<10;i++)
    {
        y[i]=x[i]+dw*(2.0*rand()/RAND_MAX-1);
        if(y[i]>mue[i]+dw)y[i]-=2*dw;else if(y[i]<mue[i]-dw)y[i]+=2*dw;
    }
}

void MC(double X[10],double T)
{   
    double y[10];
    for(int i=0;i<MCS;i++)
    {
        change(X,y);
        double de=0;double prof=0;
        prof=product(f(y))-costs(y);
        de=profit-prof;//old-new ; 因为寻找最大值
        if(de<0||exp(-de/T)*RAND_MAX>rand())
        {
            profit=prof;
            for(int j=0;j<10;j++)
            {
                X[j]=y[j];
            }
            if(profit>profit_best)
            {
                profit_best=profit;
                for(int k=0;k<10;k++)
                {
                     X_bestnow[k]=y[k];
                }
            }
        }



    }
}


int main()
{   
    srand((unsigned)time(NULL));
    initial(X);
    double T=T0;
    for(int i=0;i<100;i++)
    {
        T-=i*dt;
        MC(X,T);
    }
    for(int i=0;i<10;i++)
    {
        cout<<X_bestnow[i]<<' '<<endl;
        
    }
    cout<<profit_best<<endl;


    system("pause");
    return 0;
}