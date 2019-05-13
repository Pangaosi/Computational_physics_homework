#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<time.h>

using namespace std;

const double Pi=4*atan(1);
const double MildRand=100000000.0/RAND_MAX;//
const int ToT=10000;//总时间间隔数量
const double dt=1e-4;
const int L=1<<6;
const double k;//Boltzman 
const double dd;//晶格间距
// double A[4][L][L] 数组存放晶格内原子位置以及速度信息 
// A[0][i][j] vx, A[1][i][j] vy; A[2][i][j] x; A[3][i][j] y;


void initial(double A[4][L][L],double T0)
{
    double de=0;
    double vx,vy,x,y;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            A[0][i][j]=MildRand*rand();//速度范围从0-inf ;/RAND_MAX 是为了取得非整数值
            A[1][i][j]=MildRand*rand();
            A[2][i][j]=dd*1.0*rand()/RAND_MAX;
            A[3][i][j]=dd*1.0*rand()/RAND_MAX;
        }
    }
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            vx=MildRand*rand();//速度范围从0-inf ;/RAND_MAX 是为了取得非整数值
            vy=MildRand*rand();
            x=dd*1.0*rand()/RAND_MAX;
            y=dd*1.0*rand()/RAND_MAX;

            de=
            
        }
    }
}

void heatsides(double A[4][L][L],double Tl,double Tr,double Td,double Tu) //上下左右边温度加热至规定温度
{
    double vx,vy,x,y;
    for(int i=0;i<L;i++) //速度
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
