#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

double half(double y[],double lnyy[],int n,int N,double &n0)//y 为衰变数据，yy为传出的拟合数据，n为x数据个数，N为拟合的数据个数
{   
    double halflife;
    double lny[n];
    double k,b;//斜率和截距
    for(int j=0;j<n;j++)
    {
        lny[j]=log(y[j]); 
    }
    double kk[4]={0,0,0,0};
    for(int m=0;m<n;m++)
    {
        kk[0]+=m;
        kk[1]+=m*m;
        kk[2]+=lny[m];
        kk[3]+=m*lny[m];
    }
    b=(kk[0]*kk[3]-kk[1]*kk[2])/(kk[0]*kk[0]-(n)*kk[1]);
    k=(kk[0]*kk[2]-(n)*kk[3])/(kk[0]*kk[0]-(n)*kk[1]);
    double xx;//拟合用横坐标
    const double tt=double(n)/N;
    for(int i=0;i<N;i++)
    {   
        xx=i*tt;
        lnyy[i]=k*xx+b;
    }
    halflife=-log(2)/k;
    n0=exp(b)/(exp(b)-1);
    return halflife;
}



int main()
{
    int n;
    cout<<"data number:"<<endl;cin>>n;
    cout<<"data:"<<endl;
    double y[n];//创建数组
    double halflife;//半衰期
    double n0;//初始原子核数量
    for(int p=0;p<n;p++)
    {
        cin>>y[p];
    }
    cout<<"拟合个数"<<endl;
    int N;//拟合个数
    cin>>N;
    double yy[N];
    
    halflife=half(y,yy,n,N,n0);

    ofstream file1,file2;
    file1.open("lnylnt.dat",ios::out);
    file2.open("origindata.dat",ios::out);
    double tt;
    tt=double(n-1)/N;
    for(int k=0;k<N;k++)
    {
        file1<<1+k*tt<<'\t'<<yy[k]<<endl;  // ln[y] 插值
    }
    for(int q=0;q<n;q++)
    {
        file2<<q+1<<'\t'<<y[q]<<endl;  //原始数据,不做任何处理
    }
    file1.close();file2.close();

    cout<<"halflife:"<<halflife<<"T"<<endl;// T为时间段


    cout<<"success"<<endl;
    system("pause");
    return 0;
}