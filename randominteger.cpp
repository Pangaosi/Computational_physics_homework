// 8. random integer 0-1-2  25% 50% 25%

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<time.h>

using namespace std;

//srand((unsigned)time(NULL));

int g_012()
{
    int re=0;
    re=4*rand()/RAND_MAX;
    if(re<2)
    return re;
    else
    return re-1;
}

int main()
{
    srand((unsigned)time(NULL));
    double n0,n1,n2;// number of 0 1 2 
    double N=0;//total number= n0+n1+n2
    n0=0;n1=0;n2=0;
    int re=0;
    const int M=10000;//step number
    for(int i=0;i<M;i++)
    {   N++;
        re=g_012();
        cout<<re<<endl;
        if(re==0)
        {n0++;}
        else if(re==1)
        {n1++;}
        else
        {n2++;}
    }

    cout<<"---------------"<<endl;
    cout<<"n0 "<<"n1 "<<"n2 "<<endl;
    cout<<n0/N<<' '<<n1/N<<' '<<n2/N<<endl;




    system("pause");
    return 0;
}