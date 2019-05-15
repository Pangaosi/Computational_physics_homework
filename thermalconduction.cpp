#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<time.h>
#include<fftw3.h>

using namespace std;

//蒙卡参数
const double Pi=4*atan(1);
//const double MildRand=100000000.0/RAND_MAX;//
const int ToT=10000000;//��ʱ������� 总时间间隔数
const double dw=2000.0;
const double dt=1e-5;
const int L=1<<6;

// 晶格参数
const double m=1e-26;//质量
const double kb=1.38e-23;//Boltzman 
const double dd=1e-11;//������ 晶格距离一半
const double k2=0.5;//FPU势能参数 k/2
const double b4=0.3/4;//FPU势能参数 b/4 V（r)=k/2*r*r+b/3*x*x*x*x;
 //double A[2][L] ;//�����ž�����ԭ��λ���Լ��ٶ���Ϣ  
// A[0][i] vx, A[1][i] x; 

//分子动力学参数
const double h=1e-13;//求倒间隔

double derv(double (*f)(double),double x)
{   
    return (f(x+h)-f(x-h))/(2*h);
}

double derv2(double (*f)(double),double x)
{
    return (f(x+h)+f(x-h)-2*f(x))/(h*h);
}

double derv3(double (*f)(double),double x)
{
    
    return (f(x+h)-f(x)-h*derv(f,x)-0.5*h*h*derv2(f,x))/(h*h*h);
}


double H0(double r) //  H=V(r) 能量
{
    return k2*r*r+b4*r*r*r*r;
}

void log(double accept,double T)
{
    ofstream file1("log.dat",ios::app);
    
    file1<<"<1> the rate of acception: "<< 100.0*accept<<"% "<<"Temperature: "<<T<<endl;
    file1.close();
    
}

void Randdis(double A[2][L])        // 随即初始化
{
   for(int p=0;p<L;p++)    //每个点初始化
    {
        
        A[0][p]=40000.0*rand()/RAND_MAX-20000.0;// v
        A[1][p]=dd*(2.0*rand()/RAND_MAX-1); //晶格中心为坐标系O  , x
    }  
}

void initial(double A[2][L],double T0)  //按照波尔兹曼分布初始化,并且记录log
{
    double de=0;
    double vx,x;
    double ek=0;
    int accept=0;
    Randdis(A);  
    

    for(int j=0;j<20000;j++)  //使链按照温度T0的波尔兹曼分布
    {
        for(int q=0;q<L;q++)
        {
            int i;
            i=rand()%L;
            double v0=A[0][i];
            double rr0=A[1][(i+1)%L]-A[1][i]+2*dd;
            double rl0=A[1][i]-A[1][(i-1+L)%L]+2*dd;
            vx=A[0][i]+dw*(2.0*rand()/RAND_MAX-1);
            if(vx>20000)vx-=40000;if(vx<-20000)vx+=40000;//�ٶȷ�Χ��0-inf ;/RAND_MAX ��Ϊ��ȡ�÷�����ֵ
           
            x=A[1][i]+dd*(2.0*rand()/RAND_MAX-1);
            if(x>dd) x-=2*dd;else if(x<-dd) x+=2*dd;   
       
            double rr1=A[1][(i+1)%L]-x+2*dd;
            double rl1=x-A[1][(i-1+L)%L]+2*dd;
        
            de=0.5*m*(vx*vx-v0*v0)-H0(rr0)-H0(rl0)+H0(rr1)+H0(rl1);  
        
            if(de<0||exp(-de/kb/T0)*RAND_MAX>rand())
            {
                A[0][i]=vx;
                A[1][i]=x;
                accept++;
            }
            
        }
        if(j>4000)
            {
                for(int pp=0;pp<L;pp++)
                ek+=0.5*m*A[0][pp]*A[0][pp];
            }
    }
    log(1.0*accept/20000/L,ek*2.0/kb/16000.0/L);
    
}

void heat(double A[2][L],int i,double T1) // 加热第i个原子至 T1
{
    double vx,x;
    double v0,x0;
    v0=A[0][i];
    x0=A[1][i];
    vx=v0+dw*(2.0*rand()/RAND_MAX-1.0);
    if(vx>20000)vx-=40000;if(vx<-20000)vx+=40000;
    x=x0+dd*(2.0*rand()/RAND_MAX-1);
    if(x>dd) x-=2*dd;else if(x<-dd) x+=2*dd;  

    double rr0=A[1][(i+1)%L]-A[1][i]+2*dd;
    double rl0=A[1][i]-A[1][(i-1+L)%L]+2*dd;

    double rr1=A[1][(i+1)%L]-x+2*dd;
    double rl1=x-A[1][(i-1+L)%L]+2*dd;

    
    double de=0.5*m*(vx*vx-v0*v0)-H0(rr0)-H0(rl0)+H0(rr1)+H0(rl1);

    if(de<0||exp(-de/kb/T1)*RAND_MAX>rand())
    {
        A[0][i]=vx;
        A[1][i]=x;
    }


}


void MoleDynamics(double A[2][L],double T0,double Tl)// 首温度 T0，尾温度Tl
{
    double dt1=dt;
    double dt2=dt1*dt/2.0;
    double dt3=dt2*dt/3.0;
    double dt4=dt3*dt/4.0;
    double x0[L-2];
    for(int t=0;t<ToT;t++)
    {   
        heat(A,0,T0);
        heat(A,L,Tl);
        double x[L-2],v[L-2];
        for(int p=1;p<L-1;p++)
        {   
            v[p-1]=A[0][p];
            x[p-1]=A[1][p];
            double rr=A[1][(p+1)%L]-x[p-1]+2*dd;
            double rl=x[p-1]-A[1][(p-1+L)%L]+2*dd;
            double x1=v[p-1];
            double v1=1.0/m*(-derv(H0,rr)+derv(H0,rl));
            double v2=1.0/m*(-derv2(H0,rr)+derv2(H0,rl));
            x[p-1]=A[1][p]+x1*dt1+v1*dt2;
            v[p-1]=A[0][p]+v1*dt1+v2*dt2;
        }
        for(int q=1;q<L-1;q++)
        {
            A[0][q]=v[q-1];
            A[1][q]=x[q-1];
        }

    }
}

int main()
{
    srand((unsigned)time(NULL));
   
   double A[2][L];
   initial(A,200);
    double ek=0;
    MoleDynamics(A,300,100);
  
   

//  system("pause");
    return 0;
}
