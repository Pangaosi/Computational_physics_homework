#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

double half(double y[],double lnyy[],int n,int N,double &n0)//y Ϊ˥�����ݣ�yyΪ������������ݣ�nΪx���ݸ�����NΪ��ϵ����ݸ���
{   
    double halflife;
    double lny[n];
    double k,b;//б�ʺͽؾ�
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
    double xx;//����ú�����
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
    double y[n];//��������
    double halflife;//��˥��
    double n0;//��ʼԭ�Ӻ�����
    for(int p=0;p<n;p++)
    {
        cin>>y[p];
    }
    cout<<"��ϸ���"<<endl;
    int N;//��ϸ���
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
        file1<<1+k*tt<<'\t'<<yy[k]<<endl;  // ln[y] ��ֵ
    }
    for(int q=0;q<n;q++)
    {
        file2<<q+1<<'\t'<<y[q]<<endl;  //ԭʼ����,�����κδ���
    }
    file1.close();file2.close();

    cout<<"halflife:"<<halflife<<"T"<<endl;// TΪʱ���


    cout<<"success"<<endl;
    system("pause");
    return 0;
}