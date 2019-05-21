#include<iostream>
#include<cmath>
#include<fstream>
#include<cstdlib>
#include<ctime>

using namespace std;




//1. 仿真步长
const int ToT=20000000;
const int PoT=10000000;
const int NoT=100;

//2.材料性质
const double J=1;
const int N=6;

int PRN()
{
   // return 2*(rand()%2)-1;
    
    int i=rand()%23;
    if(i>4)
    return -1;
    else
    return 1;
}

double Energy(int A[N][N])
{   
    double e=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            e+=(A[(i+1)%N][j]+A[i][(j+1)%N])*A[i][j];
        }
    }
    e*=-J;
    return e;
}

double Magnitization(int A[N][N])
{
    double m=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            m+=A[i][j];
        }
    }
    return 1.0*m/N/N;
}

int Descend(int A[N][N])
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if((A[(i+1)%N][j]+A[(i-1+N)%N][j]+A[i][(j+1)%N]+A[i][(j-1+N)%N])*A[i][j]>0)
            {
                A[i][j]*=-1;
                return 0;
            }
        }
    }
    return 0;
}


void Log(int A[N][N],int t)
{

    ofstream file1("CellIsing.dat",ios::app);
    file1<<"<"<<t<<">---------------------energy: "<<Energy(A)<<"----------------------Magnitization: "<<Magnitization(A)<<endl;
    for( int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            file1<<A[i][j]<<' ';
        }
        file1<<endl;
    }
    file1.clear();
}


void initial(int A[N][N],int t)
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            A[i][j]=PRN();
        }
    }
    Log(A,t);
}

int Search(double F[NoT],double e)
{   
    int flag;
    int ee=int(e)+72;//e:-80~80, ee:0~1000
    if (F[ee]==1)
    {
        flag=0;
    }
    else
    {
        flag=1;
        F[ee]=1;
    }
    return flag;
    
}

void Cell(int A[N][N])
{
    int N0=N/2;
    int e0=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {   
            if ((i+j)%2==0)
            {
                e0=A[(i+1)%N][j]+A[(i-1+N)%N][j]+A[i][(j+1)%N]+A[i][(j-1+N)%N];
                if(e0==0)
                {
                     A[i][j]*=-1*(2*(rand()%2)-1);
                }
            }
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {   
            if ((i+j)%2==1)
            {
                e0=A[(i+1)%N][j]+A[(i-1+N)%N][j]+A[i][(j+1)%N]+A[i][(j-1+N)%N];
                if(e0==0)
                {
                     A[i][j]*=-1*(2*(rand()%2)-1);
                }
            }
        }
    }
}

void Run(int A[N][N])
{   
    ofstream file2("IsingData.dat",ios::out);
    
    double Ek[NoT];
    for(int i=0;i<NoT;i++)
    {
        Ek[i]=0;
    }

    for(int p=0;p<NoT;p++)
    {  double  m=0;
        int flag;
    //    Descend(A);
    initial(A,p);
    flag=Search(Ek,Energy(A));
    if(flag==1)    
    {        for(int i=0;i<ToT;i++)
            {
                Cell(A);

                if(i>PoT-1)
                {
                    m+=Magnitization(A);
                }
            }
            file2<<Energy(A)<<' '<<1.0*m/(ToT-PoT)<<endl;
    }
        
    }
    file2.close();
}



int main()
{
    srand((unsigned)time(NULL));
    
    int A[N][N];
/*    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            A[i][j]=-1;
        }
    }
*/
/*
    initial(A,0);    
    cout<<Energy(A)<<endl;
    for(int i=0;i<10;i++)
    { double m=0;
       for(int j=0;j<ToT;j++)
       {
            Cell(A);
            m+=Magnitization(A);
        }
        cout<<1.0*m/ToT<<endl;
    } */
    Run(A);
    return 0;
}