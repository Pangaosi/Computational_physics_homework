#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;
double (*func)(double x);//����ָ�� ȫ�ֱ���

double laginte3(double (*f)(double x),double x0,double x1) //�����������ղ�ֵ ����
{
    const int N=500;//��׼ȥ��ָ�
    int n;//ʵ������ָ�
    double S=0;
    double xx;//С���䳤��
    if((x1-x0)>5)  //���ݻ��ֳ��ȷָ�����
    {
        n=(x1-x0)/5*N;
    }
    xx=(x1-x0)/n;
    for(int i=0;i<n;) //����
    {
        S+=xx*(f(x0+xx*i)+4*f(x0+xx*(i+1))+f(x0+xx*(i+2)));
        i+=2;
    }
    return S/3.0;
} 

double function(double x)  //����Ϊ�κ���������
{
    double ss;
    ss=exp(x);
    return ss;
}



int main()
{
    double s,x0,x1;
    cout<<"x0 x1"<<endl;
    cin>>x0>>x1;
    s=laginte3(function,x0,x1);
    cout<<s;


    system("pause");
    return 0;
}