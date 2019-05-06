#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

double derv(double (*f)(double),double x)
{   const double h=0.0001;
    return (f(x+h)-f(x-h))/(2*h);
}

double derv2(double (*f)(double),double x)
{
    const double h=0.0001;
    return (f(x+h)+f(x-h)-2*f(x))/(h*h);
}

double derv3(double (*f)(double),double x)
{
    const double h=0.0001;
    return (f(x+h)-f(x)-h*derv(f,x)-0.5*h*h*derv2(f,x))/(h*h*h);
}


void RungeKutta(double I[],double Ii[],int N,double I0,double Ii0,double R,double L,double C,double (*f)(double x)) // didt2+R/L*didt+i/(LC)=dUdt;//该模式为
{
    const double tt=2*atan(1)/1000;
    const double dt1=tt;
    const double dt2=dt1*tt/2.0;
    const double dt3=dt2*tt/3.0;
    const double dt4=dt3*tt/4.0;
    const double a=R/L;const double b=1/(L*C);

    I[0]=I0;Ii[0]=Ii0;
    for(int i=1;i<N+1;i++)
    {   
    double t=tt*i;
    double dI1=Ii[i-1];//一阶导数,初值
    double dI2=f(t)-a*dI1-b*I[i-1];
    double dI3=derv(f,t)-a*dI2-b*dI1;
    double dI4=derv2(f,t)-a*dI3-b*dI2;
    double dI5=derv3(f,t)-a*dI4-b*dI3;
    I[i]=I[i-1]+dI1*dt1+dI2*dt2+dI3*dt3+dI4*dt4;
    Ii[i]=Ii[i-1]+dI2*dt1+dI3*dt2+dI4*dt3+dI5*dt4;
    }

    
}

double f(double x)
{
    return 0;
}

int main()
{   
    int N=1000;
    double I[N*30],Ii[N*30];
    double R,L,C;
    cout<<"R,L,C"<<endl;
    cin>>R>>L>>C;
    double I0,Ii0;
    cout<<"I,dI/dt初始值";
    cin>I0>>Ii0;
    RungeKutta(I,Ii,N*30,I0,Ii0,0,R,L,C,f);

    ofstream file1;
    file1.open("RLC.dat",ios::out);
    for(int i=0;i<30*N+1;i++)
    {
        file1<<2*atan(1)/(30*N)*i<<'\t'<<I[i]<<endl;
    }

    file1.close();
    cout<<"success";

    system("pause");
    return 0;
}