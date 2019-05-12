#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>

using namespace std;
const int MN=50000;//蒙卡总次数
const int MNP=10000;//前期蒙卡
const double dw=1;//窗口
const int L=1<<4;//L*L
double T;//temperature
const double M[3]={0,0,0};//磁化强度
const double J=1;
const double dd=2e-22;
const double hz=0.01;//z方向h
const double kb=1;
double Dr[3]={1,0,0};
double Dl[3]={-1,0,0};
double Du[3]={0,1,0};
double Dd[3]={0,-1,0};

double x1[3]={0,0,0};double x2[3]={0,0,0};double x3[3]={0,0,0};double x4[3]={0,0,0};
double x5[3]={0,0,0};double x6[3]={0,0,0};double x7[3]={0,0,0};double x8[4]={0,0,0};



const double Pi=4.0*atan(1);

class spin   //自旋 状态  类
{
    public:
    double theta,phi;  //球坐标
    spin()
    {
        theta=0;
        phi=0;
    }
};

void initial(int n,int m,spin *sp) //n*m 2 dimention
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {
            //int k=i*m+j;
            (sp+i*m+j)->phi=2*Pi*rand()/RAND_MAX;  //对 phi 随即初始化
            sp[i*m+j].theta=acos(2.0*rand()/RAND_MAX-1); // 对 cos(theta) 随即初始化
        }
    }
}

void rpt2xyz(spin sp,double *xyz)// 球坐标转笛卡尔坐标
{
    double theta,phi;
    theta=sp.theta;phi=sp.phi;
    xyz[0]=1.0*sin(theta)*cos(phi);
    xyz[1]=1.0*sin(theta)*sin(phi);
    xyz[2]=1.0*cos(theta);
}

double dot(spin sp1,spin sp2) //点乘
{
    double xyz1[3]={0,0,0};
    double xyz2[3]={0,0,0};
    rpt2xyz(sp1,xyz1);
    rpt2xyz(sp2,xyz2);
    return xyz1[0]*xyz2[0]+xyz1[1]*xyz2[1]+xyz1[2]*xyz2[2];
}

double dot(double *xyz1,double *xyz2) // 重定义
{
    return xyz1[0]*xyz2[0]+xyz1[1]*xyz2[1]+xyz1[2]*xyz2[2];
}

double *cross(spin sp1,spin sp2,double *xyz3)
{
    double xyz1[3]={0,0,0};
    double xyz2[3]={0,0,0};
   // double xyz3[3]={0,0,0};
    rpt2xyz(sp1,xyz1);
    rpt2xyz(sp2,xyz2);
    double x1=xyz1[0];double x2=xyz1[1];double x3=xyz1[2];
    double y1=xyz2[0];double y2=xyz1[1];double y3=xyz1[2];
    xyz3[0]=x2*y3-x3*y2;
    xyz3[1]=-x1*y3+x3*y1;
    xyz3[2]=x1*y2-x2*y1;
    return xyz3;
}

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

double scpro(double *D,spin sp1,spin sp2) //标量三重积
{
 double sum=0;
 rpt2xyz(sp1,x1);
 rpt2xyz(sp2,x2);
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            for(int k=0;k<3;k++)
            {
                sum+=levi(i+1,j+1,k+1)*D[i]*x1[j]*x2[k];
            }
        }
    }
    return sum;
}


/*
double *D(int i,int j,int p,int q,double *xyz)//(i,j,0)->(p,q,0)
{
    double *Dij=xyz;
    
    Dij[0]=p-i;
    Dij[1]=q-j;
    Dij[2]=0;
    return Dij;
}
*/
void Magnetization(spin *sp,double *Mxyz)
{   
    Mxyz[0]=0;Mxyz[1]=0;Mxyz[2]=0;
    double xyz[3]={0,0,0};
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {
            int k=i*L+j;
            rpt2xyz(sp[k],xyz);
            Mxyz[0]+=xyz[0];
            Mxyz[1]+=xyz[1];
            Mxyz[2]+=xyz[2];   
        }
    }
    Mxyz[0]/=L*L;
    Mxyz[1]/=L*L;
    Mxyz[2]/=L*L;
}
//double Dd[3]={0,-1,0};
double energy(spin *sp)
{
    double sum=0;
    for(int i=0;i<L;i++)
    {
        for(int j=0;j<L;j++)
        {   
            int k=i*L+j;
            double xyz[3]={0,0,0};double x1[3]={0,0,0};double x2[3]={0,0,0};
            double x4[3]={0,0,0};double x5[3]={0,0,0};
           
            sum+=J*(dot(sp[k],sp[L*i+(j+1)%L])+dot(sp[k],sp[(k+L)%(L*L)]));//J*(Si.Sj)
            sum+=dd*(cross(sp[k],sp[L*i+(j+1)%L],x4)[0]+cross(sp[k],sp[(k+L)%(L*L)],x5)[1]);
        //    sum+=dd*(dot(Dr,cross(sp[k],sp[L*i+(j+1)%L],x4))+dot(Du,cross(sp[k],sp[(k+L)%(L*L)],x5)));
            rpt2xyz(sp[k],xyz);
            sum+=hz*xyz[2];
        }
    }
    return sum/L/L;
}
/*
void save(int n,int m,const spin *sp)//n*m,xyz
{
    ofstream file1,file2,file3;
    file1.open("spinx.dat",ios::app);
    double xyz[3];
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
        {   
            int k=m*i+j;
            rpt2xyz(sp[k],xyz);       
            file1<<j<<' '<<i<<' '<<xyz[0]<<' '<<xyz[1]<<' '<<xyz[2]<<endl;
        }
    }
    file1.close();
}
*/
void MC(spin* sp)
{   
    double Me[3]={0};
    double Se=0;double Se2=0;//Me 磁化强度，Se总能量
    double accept=0;//接受率
    for(int p=0;p<MN;p++)
    {
        for(int q=0;q<L*L;q++)
        { 
            int k=rand()%(L*L); 
            int i=k/L;int j=k%L;//随即取点
            spin sp0;//临时变量sp0
            sp0.phi=sp[k].phi+Pi*dw*(2.0*rand()/RAND_MAX-1);
            double costheta;
            costheta=cos(sp[k].theta)+dw*(2.0*rand()/RAND_MAX-1);
            if(costheta>1) costheta=2-costheta;
            if(costheta<-1) costheta=-2-costheta;
            sp0.theta=acos(costheta); //生成下一采样   
            
          //  
            
            double de=0;//能量变化 new-old
            de-=J*(dot(sp[k],sp[L*i+(j+1)%L])+dot(sp[k],sp[(k+L)%(L*L)])+dot(sp[k],sp[L*i+(L+j-1)%L])+dot(sp[k],sp[(k-L+L*L)%(L*L)]));
            de-=cross(sp[k],sp[L*i+(j+1)%L],x5)[0]+cross(sp[k],sp[(k+L)%(L*L)],x6)[1]-cross(sp[k],sp[L*i+(L+j-1)%L],x7)[0]-cross(sp[k],sp[(k-L+L*L)%(L*L)],x8)[1];
            de-=hz*cos(sp[k].theta);
            de+=J*(dot(sp0,sp[L*i+(j+1)%L])+dot(sp0,sp[(k+L)%(L*L)])+dot(sp0,sp[L*i+(L+j-1)%L])+dot(sp0,sp[(k-L+L*L)%(L*L)]));
            de+=cross(sp0,sp[L*i+(j+1)%L],x5)[0]+cross(sp0,sp[(k+L)%(L*L)],x6)[1]-cross(sp0,sp[L*i+(L+j-1)%L],x7)[0]-cross(sp0,sp[(k-L+L*L)%(L*L)],x8)[1];
            de+=hz*cos(sp0.theta);

            if(de<0||exp(-de/(kb*T))*RAND_MAX>rand())
            {
                sp[k].phi=sp0.phi;
                sp[k].theta=sp0.theta;
                accept++;
            }
        }

        if(p>MNP-1)
        {
            double m[3];
            Magnetization(sp,m);
            Me[0]+=m[0];
            Me[1]+=m[1];
            Me[2]+=m[2];
            Se+=energy(sp);
            Se2+=energy(sp)*energy(sp);
        }
    }
    Me[0]/=MNP;Me[1]/=MNP;Me[2]/=MNP;
    Se/=MNP;
    Se2/=MNP;
    double cp=(Se2-Se*Se)/(L*L*kb*T*T);
    ofstream file1;
    file1.open("Heisenberg.dat",ios::app);
    file1<<T<<' '<<hz<<' '<<Se<<' '<<cp<<' '<<Me[0]<<' '<<Me[1]<<' '<<Me[2]<<' '<<1.0*accept/MN/(L*L)<<endl;
    cout<<T<<' '<<hz<<' '<<Se<<' '<<cp<<' '<<Me[0]<<' '<<Me[1]<<' '<<Me[2]<<' '<<1.0*accept/MN/(L*L)<<endl;
    file1.close();
}

int main()
{
    srand((unsigned)time(NULL));
    int LL=L*L;
    spin sp[LL];
    initial(L,L,sp);

    for(T=10;T>2;T-=0.5)
    {
        cout<<T<<endl;
        MC(sp);
    }

/*
    double x[3]={0,0,0};
    ofstream file1;
    file1.open("text.dat",ios::out);
    for(int i=0;i<10000;i++)
    {
        sp.phi=2*Pi*rand()/RAND_MAX;  //对 phi 随即初始化
        sp.theta=acos(2.0*rand()/RAND_MAX-1);
        rpt2xyz(sp,x);
        file1<<x[0]<<' '<<x[1]<<' '<<x[2]<<endl;
    } 
    double x1[3]={0,1,1};
    double x2[3]={1,0,-2};
    spin spp;
    spp.phi=Pi/3;
    spp.theta=Pi/6;
    rpt2xyz(spp,x);
    cout<<x[0]<<' '<<x[1]<<' '<<x[2]<<endl;
    cout<<dot(x1,x2)<<endl; */
 //   file1.close();
    return 0;
}