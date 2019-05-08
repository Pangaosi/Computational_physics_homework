#include<iostream>
#include<fstream>
#include<cmath>
#include<fftw3.h>
#include<complex>

using namespace std;

const double Pi=4*atan(1);
const int m=7;
const int N=1<<m;


int main()
{
    fftw_complex *in;
    fftw_complex *out;
    fftw_plan p;
    in=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
    out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);

    for(int i=0;i<N;i++)
    {  double x=i*2.0*Pi/N;
        in[i][0]=sin(x)+sin(25*x+1)+sin(52*x+2);
        in[i][1]=0;
    }
 /*   for(int i=0;i<N;i++)  //initialize
    {
        for(int j=0;j<N;j++)
        {   
            int k=i*N+j;
            double x,y;
            x=2*Pi/double(N)*double(i);
            y=2*Pi/double(N)*double(j);
            in[k][0]=sin(5*x+7*y);
            in[k][1]=0;
     
        }
    } 
*/

    p=fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);

 //   ofstream OUTRE,OUTIM;
 //   OUTRE.open("fftre.dat",ios::out);
 //   OUTIM.open("fftim.dat",ios::out);


 //   const int L=(N-1-(N-1)%2)/2;
 /*   for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
           int k=N*i+j;
            OUTRE<<2*Pi/double(N)*double(j)<<' '<<2*Pi/double(N)*double(i)<<' '<<out[k][0]<<endl;
            OUTIM<<2*Pi/double(N)*double(j)<<' '<<2*Pi/double(N)*double(i)<<' '<<out[k][1]<<endl;
            
            
        }
    
        OUTRE<<endl;
        OUTIM<<endl;
    }
    OUTRE.close();
    OUTIM.close(); 
    */
   
    ofstream file1;
    file1.open("DFT-fft.dat",ios::out);
    for(int k=0;k<N;k++)
    {
        file1<<k<<' '<<out[k][0]<<' '<<out[k][1]<<endl;
    }
    file1.close();
    fftw_free(in);
    fftw_free(out);


    
    return 0;
}