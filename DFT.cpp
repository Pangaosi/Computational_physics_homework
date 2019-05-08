#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

const double Pi= 4*atan(1);

class complex
{
    public:
    double re;
    double im;
    complex()
    {
        re=0;
        im=0;
    }
};

void DFT(int N,complex *F,complex *W)
{
    const double aa=1.0/sqrt(N);   
    for(int j=0;j<N;j++)
    {   
        const double co=2*Pi*j/double(N);//coeffient
        for(int k=0;k<N;k++)
        {
           double re=(F+k)->re;
           double im=(F+k)->im;
           (W+j)->re+=(cos(co*k)*re+sin(co*k)*im);
           (W+j)->im+=(cos(co*k)*im-sin(co*k)*re);
        }
        (W+j)->re*=aa;
        (W+j)->im*=aa;
    }     
} 

void   trans(int n,int m,complex *W,complex *Wt)//ԭ��Ϊn*m ת��Ϊ m*n
{
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            //Wt[i][j]  i*n+j W[j][i]  j*m+i
            Wt[i*n+j].re=W[j*m+i].re;
            Wt[i*n+j].im=W[j*m+i].im;
        }
    }

}

void showRE(int n,int m,complex *W) //��ӡ RE
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            cout<<W[i*m+j].re<<' ';
        }
        cout<<endl;
    }
}


void DFT_2D(int N,int M,complex *F,complex *W)
{
    for(int i=0;i<N;i++)
    {
        DFT(M,F+M*i,W+M*i);     
    }
    int NM=N*M;
    complex Wt[NM];
    trans(N,M,W,Wt);
    complex Wtt[NM];
    for(int j=0;j<M;j++)
    {
        DFT(N,Wt+N*j,Wtt+N*j);
    }
    trans(M,N,Wtt,W);
}

int main()
{   
    int m=7;
    int N=1<<m;
    int NN=N*N;
/*    complex F[N];
    for(int i=0;i<N;i++)
    {  double x=i*2.0*Pi/N;
        F[i].re=sin(x)+sin(25*x+1)+sin(52*x+2);
        F[i].im=0;
    } 
    */
    complex F2[NN];complex W2[NN];
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            double x=2*Pi/double(N)*double(j);
            double y=2*Pi/double(N)*double(i);
            F2[N*i+j].re=sin(5.0*x+7.0*y);
            F2[N*i+j].im=0;
        }
    } 
    DFT_2D(N,N,F2,W2);
    

    ofstream file1;
    file1.open("DFT_22D.dat",ios::out);
    for(int k=0;k<N;k++)
    {
        for(int p=0;p<N;p++)
        {   
            double x=2*Pi/double(N)*double(p);
            double y=2*Pi/double(N)*double(k);
            file1<<x<<' '<<y<<' '<<W2[k*N+p].re<<' '<<W2[k*N+p].im<<' '<<F2[k*N+p].re<<endl;
        }
        
    }
    file1.close();
    
/*
    complex W[N];
    DFT(F,N,W);
    ofstream file1;
    file1.open("DFT.dat",ios::out);
    for(int k=0;k<N;k++)
    {
        file1<<k<<' '<<W[k].re<<' '<<W[k].im<<endl;
    }
    file1.close();
*/   

 
    
 return 0;   
}
