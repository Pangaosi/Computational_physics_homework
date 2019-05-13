#include<iostream>
#include<fstream>

using namespace std;

extern "C"
{
    extern int dsyev_(char*,char*,int*,double*,int *,double*,double*,int*,int*);
}

int main()
{
    double data[100];
    ifstream file1;
    file1.open("matrix.dat",ios::in);

    for(int i=0;i<10;i++)
    {
        for(int j=0;j<10;j++)
        {
            file1>>data[10*i+j];
        }
    }
    file1.close();

    char JOBZ='V';
    char UPLO='L';
    int n=10;
    int LDA=n;
    double *W=new double[n];
    int lwork=6*n;

    double *work=new double[lwork];

    int info;

    dsyev_(&JOBZ,&UPLO,&n,data,&LDA,W,work,&lwork,&info);

    cout<<"info: "<<info<<endl;
    cout<<"------eigenvalues--------"<<endl;
    for(int i=0;i<10;i++)
    {
        cout<<W[i]<<' ';
    }
    cout<<endl;

    cout<<"-------eigenvectors--------"<<endl;
    for(int i=0;i<10;i++)
    {
        cout <<"<"<<i<<">:";
        for(int j=0;j<10;j++)
        {
            
            cout<<data[j*n+i]<<' ';
        }
        cout<<endl;
    }
    delete []W;
    delete []work;
    return 0;
}
