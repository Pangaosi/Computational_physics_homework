#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

int levi(int i,int j,int k)
{
    if ((j-i)*(k-j)*(i-k)==0)
    {
        return 0;
    }
    else if ((j-i)*(k-j)*(i-k)<0)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}


double det(double A[3][3])//3行列式
{   
    double sum=0;
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                sum+=levi(i+1,j+1,k+1)*A[0][i]*A[1][j]*A[2][k];
            }
        }
    }
    return sum;
}

void copy(double A[3][3],double B[3][3])
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            B[i][j]=A[i][j];
        }
    }
}

void showtime(double A[3][3])
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            cout<<A[i][j]<<' ';
        }
        cout<<endl;
    }
}

void solve(double A[3][3],double *B,double *x) //3阶求解
{   double AA[3][3];
    copy(A,AA);
    if (det(AA)!=0)
    {
        for(int i=0;i<3;i++)// 从行开始高斯消元
        {   double a=AA[i][0];
            for(int j=0;j<3;j++)
            {
                AA[i][j]/=a;
            }
            B[i]/=a;
        }
        for(int p=0;p<3-1;p++)
        {
            for(int q=0;q<3;q++)
            {
                AA[p+1][q]-=AA[0][q];
            }
            B[p+1]-=B[0];
        }
        for(int e=0;e<3-1;e++)/* code */
        {   double b=AA[e+1][1];
            for(int r=0;r<3-1;r++)
            {
                AA[e+1][r+1]/=b;
            }
            B[e+1]/=b;
        }
        for(int t=0;t<3-1;t++)
        {
            AA[2][t+1]-=AA[1][t+1];
        }
        B[2]-=B[1];

        x[2]=B[2]/AA[2][2];
        x[1]=(B[1]-x[2]*AA[1][2])/AA[1][1];
        x[0]=(B[0]-x[2]*AA[0][2]-x[1]*AA[0][1])/AA[0][0];
    }
    else
    {
        cout<<"det=0";
    }
}




int main()
{   
    double R1,R2,R3,Ra,Rs,Rx,U0;
    cout<<"R1,R2,R3,Ra,Rs,Rx,U0"<<endl;
    cin>>R1>>R2>>R3>>Ra>>Rs>>Rx>>U0;
    double R[3][3]={{Rs,R1,R2},{-Rx,R1+Rx+Ra,-Ra},{-R3,-Ra,R2+R3+Ra}};
    double V[3]={U0,0,0};
    double x[3]={0,0,0};
    solve(R,V,x);
 
    cout<<x[0]<<' '<<x[1]<<' '<<x[2]<<endl;
    
    return 0;
}