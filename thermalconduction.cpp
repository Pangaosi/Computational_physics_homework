#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<time.h>

using namespace std;

const double Pi=4*atan(1);
const double MildRand=100000000.0/RAND_MAX;//
const int ToT=10000;//��ʱ��������
const double dt=1e-4;
const int L=1<<6;
const double k;//Boltzman 
const double dd;//������
// double A[4][L][L] �����ž�����ԭ��λ���Լ��ٶ���Ϣ 
// A[0][i][j] vx, A[1][i][j] vy; A[2][i][j] x; A[3][i][j] y;


void initial(double A[4][L][L],double T0)
{
    double de=0;
    double vx,vy,x,y;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            A[0][i][j]=MildRand*rand();//�ٶȷ�Χ��0-inf ;/RAND_MAX ��Ϊ��ȡ�÷�����ֵ
            A[1][i][j]=MildRand*rand();
            A[2][i][j]=dd*1.0*rand()/RAND_MAX;
            A[3][i][j]=dd*1.0*rand()/RAND_MAX;
        }
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            vx=MildRand*rand();//�ٶȷ�Χ��0-inf ;/RAND_MAX ��Ϊ��ȡ�÷�����ֵ
            vy=MildRand*rand();
            x=dd*1.0*rand()/RAND_MAX;
            y=dd*1.0*rand()/RAND_MAX;

            de=
            
        }
    }
}

void heatsides(double A[4][L][L],double Tl,double Tr,double Td,double Tu) //�������ұ��¶ȼ������涨�¶�
{
    double vx,vy,x,y;
    for(int i=0;i<L;i++) //�ٶ�
    {
        vx=10000000.0*rand()/RAND_MAX;
        vy=10000000.0*rand()/RAND_MAX;
        x=dd*1.0*rand()/RAND_MAX;
        y=dd*1.0*rand()/RAND_MAX;
    }
}

int main()
{
   

    
  system("pause");
    return 0;
}
