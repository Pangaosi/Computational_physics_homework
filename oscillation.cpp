#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

int main()
{  
    const int N=5000;//��׼����
    double x0,v0,k,b,m;//��ʼλ�ã���ʼ�ٶȣ�����ϵ��������ϵ��,����
    double t;//�ܹ۲�ʱ��
    int n;
    cout<<"��ʼλ�ã���ʼ�ٶȣ�����ϵ��������ϵ��������,�۲�ʱ��[SI]"<<endl;
    cin>>x0>>v0>>k>>b>>m>>t;
    if(t>3)
        n=(t/3)*N;//ʵ�ʴ���    
    else
        n=N; 
    double tt=t/n;
    double x[n],v[n];//a[n];
    v[0]=v0;x[0]=x0;
    for(int i=0;i<n-1;i++)
    {
        v[i+1]=v[i]+(-(k/m)*x[i]-(b/m)*v[i])*tt;
        x[i+1]=x[i]+0.5*(v[i]+v[i+1])*tt;
        
    }

    ofstream file1;
    file1.open("positon.dat",ios::out);
    for(int p=0;p<n;p++)
    {   double ttt=tt*p;
        file1<<ttt<<'\t'<<x[p]<<endl;
    }
    file1.close();
    cout<<"success";
    
    system("pause");
    return 0; 
}