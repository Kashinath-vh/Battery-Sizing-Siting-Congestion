#include<iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cstdlib>
using namespace std;

class gen{
public:
int isDG,qlim,type;
double mp,nq,V0,w0,Qmax,a,b,c,Pmax,Smax,Pg,Qg,Pmin,Qper;
gen();
};

gen::gen ()
{
isDG=0;
qlim=0;
V0=1.01;
w0=1.0;
Qmax=10000000;
Smax=10000000;
type=6;
Pmin=0;
Qper=0.01;
mp=0.2;
nq=0.87;
}

class load{
public:
int type;
double Pl0,Ql0,alpha,beta,kpf,kqf;
load();
};

load::load()
{
type=0;
Pl0=0;
Ql0=0;
alpha=0;
kpf=0;
kqf=0;
beta=0;
}

class line{
public:
int isON,from,to;
double r,x,limit;
complex<double> flow;
line();
};

line::line()
{
isON=0;
limit=500;
}

double w0,PI=3.141592654,w=1.0,V0=1.0,Ploss,tolerance,base_V,base_VA,base_I,pf,loaddatai[101],loaddatar[101],J[250][250],J1[250][250],
battpower[101],del,Vsort[120][3],Vmin=0.95,Vmax=1.05;
int no_lines,no_tie,no_iterations,no_buses,location,LBN,no_gens,no_DG,m,n,found[120];

gen DG[120];
load LD[120];
line tran[140];

complex<double> I[120],Scal[120],V[120],Z_base,SLtotal,SL[120];

int sortdescend(complex<double>* x, int n)
{
int i,j;
double temp;
for(i=1;i<=n;i++)
{
Vsort[i][1]=abs(x[i]);
Vsort[i][2]=i;
}
for(i=1;i<=n;i++)
for(j=1;j<=(n-i+1);j++)
if(Vsort[j][1]<Vsort[j+1][1])
{
temp=Vsort[j][1];
Vsort[j][1]=Vsort[j+1][1];
Vsort[j+1][1]=temp;
temp=Vsort[j][2];
Vsort[j][2]=Vsort[j+1][2];
Vsort[j+1][2]=temp;
}
return 1;
}

int findn( int element)
{
int row,i,n=0;
for(i=1;i<=(no_lines+no_tie);i++)
{
if((tran[i].from==element)&&(tran[i].isON==1))
{
found[n+1]=tran[i].to;
n++;
}
}
found[0]=n;
return n;
}

double max(double* a, int n)
{
int i;
double x;
x=a[1];
for(i=2;i<=n;i++)
if(a[i]>x)
x=a[i];
return x;
}

double min(complex<double>* a, int n)
{
int i;
double x;
x=abs(a[1]);
for(i=2;i<=n;i++)
if(abs(a[i])<x)
x=abs(a[i]);
return x;
}

double i2soc(int i)
{
double x=0;
x=(1.0/del +50.0-double(i))/(1.0/del);
return(x);
}

int soc2i(double x)
{
return(1.0/del +50.0 -round((1.0/del)*x));
}

complex<double> sum(complex<double> * x,int n)
{
int i;
complex<double> y;
for(i=1;i<=n;i++)
y=y + x[i];
return y;
}

double absmax(double* x,int n)
{int i;
double y;
y=abs(x[1]);
for(i=2;i<=n;i++)
if(abs(x[i])>y)
y=abs(x[i]);
return y;
}

int invshipley(double a[250][250],int n )
{
int i,j,k,precision_no=6;
for(k=1;k<=n;k++)
{
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
{
if(i==k || j==k)
continue;
a[i][j]=a[i][j] - (a[i][k]*a[k][j])/a[k][k];
}
for(i=1;i<=n;i++)
{if (i==k)
continue;
a[i][k]=(-1*a[i][k])/a[k][k];
}
for(j=1;j<=n;j++)
{
if (j==k)
continue;
a[k][j]=(-1.0*a[k][j])/a[k][k];
}
a[k][k]=(-1.0/a[k][k]);
}
for(i=1;i<=n;i++)
for(j=1;j<=n;j++)
a[i][j]=-1.0*a[i][j];
return 0;
}

int islandnrlf(int m, double Pdg1=0)
{
int i,j,k,l,from_bus,to_bus;
double delP[250],delV[250],Vang,Vmag,r,x;
complex<double> SG[120],z,Y[120][120];
w=1.0;
V[1]=1.0 + 0i;
for(i=2;i<=no_buses;i++)
V[i]= 1. + 0i;
for(i=1;i<=no_buses;i++)
SG[i]=complex<double>(0,0);

for(k=1;k<=no_iterations;k++)
{
//Ybus
for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
Y[i][j]=complex<double>(0,0);

for(i=1;i<=(no_lines+no_tie);i++)
if(tran[i].isON==1)
{
from_bus=tran[i].from;
to_bus=tran[i].to;
r=tran[i].r;
x=tran[i].x;
z=complex<double>(r,w*x);
Y[from_bus][to_bus]=Y[from_bus][to_bus] - Z_base/z;
Y[from_bus][from_bus]=Y[from_bus][from_bus] + Z_base/z;
Y[to_bus][to_bus]=Y[to_bus][to_bus] + Z_base/z;
Y[to_bus][from_bus]=Y[from_bus][to_bus];
}

for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
if(DG[i].type==4)
{
r=DG[i].Pg;
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==5)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=DG[i].Pg;
}
if(DG[i].type==6)
{
x=DG[i].Qg;
r=DG[i].Pg;
}
DG[i].Qmax=sqrt(DG[i].Smax*DG[i].Smax - DG[i].Pg*DG[i].Pg);
DG[i].qlim=0;
if(x>=DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}
if(x<=-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}
DG[i].Qg=x;
SG[i]=complex<double>(r,x);
}
for(i=1;i<=no_buses;i++)
{
r=LD[i].Pl0*(pow(abs(V[i]),LD[i].alpha))*(1.0+LD[i].kpf*(w-1.0));
x=LD[i].Ql0*(pow(abs(V[i]),LD[i].beta))*(1.0+LD[i].kqf*(w-1.0));
if(i<23)
SL[i]=complex<double>(r,x)*loaddatar[m];
else
SL[i]=complex<double>(r,x)*loaddatai[m];
}

SL[location]+=complex<double>(Pdg1,0);

for(i=1;i<=no_buses;i++)
for(j=1;j<=no_buses;j++)
J[i][j]=0;

for(i=1;i<=no_buses;i++)
{
complex<double> temp=0;
for(j=1;j<=no_buses;j++)
temp+=Y[i][j]*V[j];
I[i]=temp;
}

for(i=1;i<=no_buses;i++)
Scal[i]=(V[i]*conj(I[i]));

//    J11
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
{
if(i==j)
J[i-1][j-1]=-(imag(Scal[i]) + abs(V[i])*abs(V[i])*imag(Y[i][i]));
else
J[i-1][j-1]=-(abs(V[i])*abs(V[j])*abs(Y[i][j])*sin( arg(Y[i][j]) + arg(V[j]) - arg(V[i])));
}

//J21
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
J[i-2+no_buses][j-1]= real(Scal[i]) - abs(V[i])*abs(V[i])*real(Y[i][i]);
else
J[i-2+no_buses][j-1]= -abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J12
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-1][j-2+no_buses]= real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[i-1][j-2+no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[i-1][j-2+no_buses]+=0.5*abs(V[i])/DG[i].nq;
}
}
else
J[i-1][j-2+no_buses]=abs(V[i])*abs(V[j])*abs(Y[i][j])*cos(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J22
for(i=2;i<=no_buses;i++)
for(j=2;j<=no_buses;j++)
if(i==j)
{
J[i-2+no_buses][j-2+no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]) + LD[i].beta*imag(SL[i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[i-2+no_buses][j-2+no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[i-2+no_buses][j-2+no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
if(DG[i].type==4)
J[i-2+no_buses][j-2+no_buses]+=(1/DG[i].nq)*abs(V[i]);
}
}
else
J[i-2+no_buses][j-2+no_buses]=-abs(V[i])*abs(V[j])*abs(Y[i][j])*sin(arg(Y[i][j]) + arg(V[j]) - arg(V[i]));

//J13 & J23
j=2*no_buses - 1;
for(i=2;i<=no_buses;i++)
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[i-1][j]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf*w;
J[i-2+no_buses][j]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf*w;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[i-1][j]+=(1/DG[i].mp);
if((DG[i].type==2)&&(DG[i].qlim==0))
J[i-2+no_buses][j]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[i-1][j]+=(0.5/DG[i].mp);
J[i-2+no_buses][j]-=(0.5/DG[i].mp);
}
if((DG[i].type==5)&&(DG[i].qlim==0))
J[i-2+no_buses][j]-=(1/DG[i].mp);
}
}

//J14
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-1][j]=abs(V[i])*abs(V[1])*abs(Y[i][1])*cos(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));
}

//J24
j=2*no_buses;
for(i=2;i<=no_buses;i++)
{
J[i-2+no_buses][j]=-abs(V[i])*abs(V[1])*abs(Y[i][1])*sin(arg(Y[i][1]) + arg(V[1]) - arg(V[i]));
}

//J31
for(j=2;j<=no_buses;j++)
J[2*no_buses - 1][j-1]=-(abs(V[1])*abs(V[j])*abs(Y[1][j])*sin( arg(Y[1][j]) + arg(V[j]) - arg(V[1])));

//J32
for(j=2;j<=no_buses;j++)
{
J[2*no_buses - 1][j-2+no_buses]=abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J41
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-1]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*cos(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J42
i=2*no_buses;
for(j=2;j<=no_buses;j++)
{
J[i][j-2+no_buses]=-abs(V[1])*abs(V[j])*abs(Y[1][j])*sin(arg(Y[1][j]) + arg(V[j]) - arg(V[1]));;
}

//J33 & J43
i=1;
{
complex<double> temp=0+0i;
for(l=1;l<=no_buses;l++)
{
if(abs(Y[i][l]))
{
r=real(complex<double>(-1.0,0)/Y[i][l]);
temp+=(V[i]-V[l])*(Y[i][l]/w)*(complex<double>(1,0) + Y[i][l]*complex<double>(r,0))*conj(V[i]);
}
}
complex<double> Yii=0+0i;
int z;
for(z=1;z<=no_buses;z++)
Yii+=Y[i][z];
if(abs(Yii))
{
r=real(complex<double>(1.0,0)/Yii);
temp-=V[i]*(Yii/w)*(complex<double>(1,0)+complex<double>(r,0)*Yii)*conj(V[i]);
}
J[2*no_buses - 1][2*no_buses - 1]=real(temp) + LD[i].Pl0*pow(abs(V[i]),LD[i].alpha)*LD[i].kpf*w;
J[2*no_buses][2*no_buses - 1]=-imag(temp) + LD[i].Ql0*pow(abs(V[i]),LD[i].beta)*LD[i].kqf*w;
if(DG[i].isDG)
{
if(DG[i].type==1)
J[2*no_buses - 1][2*no_buses - 1]+=(1/DG[i].mp);
if((DG[i].type==2)&&(DG[i].qlim==0))
J[2*no_buses][2*no_buses - 1]-=(1/DG[i].mp);
if(DG[i].type==3)
{
J[2*no_buses - 1][2*no_buses - 1]+=(0.5/DG[i].mp);
J[2*no_buses][2*no_buses - 1]-=(0.5/DG[i].mp);
}
if((DG[i].type==5)&&(DG[i].qlim==0))
J[2*no_buses][2*no_buses - 1]-=(1/DG[i].mp);
}
}

//J34
i=1;
J[2*no_buses - 1][2*no_buses]=real(Scal[i]) + abs(V[i])*abs(V[i])*real(Y[i][i]) + real(SL[i])*LD[i].alpha;
if(DG[i].isDG)
{
if(DG[i].type==2)
J[2*no_buses - 1][2*no_buses]+=abs(V[i])/DG[i].nq;
if(DG[i].type==3)
J[2*no_buses - 1][2*no_buses]+=0.5*abs(V[i])/DG[i].nq;
}

//J44
i=1;
j=1;
J[2*no_buses][2*no_buses]= imag(Scal[i]) - abs(V[i])*abs(V[i])*imag(Y[i][i]);
if(DG[i].isDG&&(DG[i].qlim==0))
{
if(DG[i].type==1)
J[2*no_buses][2*no_buses]+=(1/DG[i].nq)*abs(V[i]);
if(DG[i].type==3)
J[2*no_buses][2*no_buses]+=(0.5/DG[i].nq)*abs(V[i]);
if(DG[i].type==4)
J[2*no_buses][2*no_buses]+=(1/DG[i].nq)*abs(V[i]);
}

for(i=1;i<=2*no_buses;i++)
for(j=1;j<=2*no_buses;j++)
J1[i][j]=J[i][j];

invshipley(J,2*no_buses);

for(i=2;i<=no_buses;i++)
delP[i-1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
for(i=no_buses;i<=2*(no_buses - 1);i++)
delP[i]=-(imag(Scal[i-no_buses+2]) + imag(SL[i-no_buses+2]) - imag(SG[i-no_buses+2]));

i=1;
delP[2*no_buses - 1]=-(real(Scal[i]) + real(SL[i]) -real(SG[i]));
delP[2*no_buses]=-(imag(Scal[i]) + imag(SL[i]) - imag(SG[i]));

for(i=1;i<=2*(no_buses);i++)
{
double temp=0;
for(j=1;j<=2*(no_buses);j++)
temp+=J[i][j]*delP[j];
delV[i]=temp;
}

for(i=1;i<=(no_buses - 1);i++)
{
Vang=arg(V[i+1])+delV[i];
Vmag=abs(V[i+1])*(1.0 + delV[i+ no_buses - 1]);
V[i+1]=polar(Vmag,Vang);
}

w+=delV[2*no_buses - 1];
Vmag=abs(V[1])*(1.0 + delV[2*no_buses]);
V[1]=polar(Vmag,0.0);

if(tolerance>absmax(delP,(2*no_buses)))  
break;
}

if(k<=no_iterations)
{
SLtotal=sum(Scal,no_buses);
Ploss=real(SLtotal)*base_VA/1000;
for(i=1;i<=no_buses;i++)
if(DG[i].isDG)
{
if(DG[i].type==1)
{
r=(1/DG[i].mp)*(DG[i].w0-w);
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==2)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==3)
{
double a,b;
a=(1/DG[i].mp)*(DG[i].w0-w);
b=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
r=(a+b)/2.0;
x=(b-a)/2.0;
}
if(DG[i].type==4)
{
r=DG[i].Pg;
x=(1/DG[i].nq)*(DG[i].V0-abs(V[i]));
}
if(DG[i].type==5)
{
x=(1/DG[i].mp)*(-DG[i].w0+w);
r=DG[i].Pg;
}
if(DG[i].type==6)
{
x=DG[i].Qg;
r=DG[i].Pg;
}
DG[i].qlim=0;
if(x>=DG[i].Qmax)
{
x=DG[i].Qmax;
DG[i].qlim=1;
}

if(x<=-DG[i].Qmax)
{
x=-DG[i].Qmax;
DG[i].qlim=1;
}
DG[i].Pg=r;
DG[i].Qg=x;
}

return 1;
}
else
cout<<"\nLoad-flow solution did not converge\n";
return 0;
}

int main()
{
int i,j,k,kmin,kmax,count,l,n1,no_hrs=24,ndel,temp,nbatt=0,LBN,e,no_to,bibc[120][120];
double r,x,price[101],cost[320][101],socmin[101],socmin1[101],socmin2[101],socmax[101],socmax1[101],socmax2[101],socst,smin,smax,temp1,soc,
socT[102],charge[102],Pd,lambda,Pdemand,ncharge=1.0,ndis=1.0,Pmax,Emax,Cbatt,cost1,Jlbn[250],
Qlbn[250],Qpenalty[250],penalty[240],tcostold,delbat,a[200],tstep,ntstep,tempsum1,tempsum2,lstep,Cq,qexit,muV[120],muIb[120],
LMP[120][27],LMPloss[120][27],LMPV[120][27],LMPI[120][27],LMPE[120][27]; 
complex<double> z,line_z[120],tempsum;
string comment;
fstream fp;
fp.open("islandip33.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
fp>>no_buses;
fp>>no_lines;
fp>>tolerance;
fp>>no_iterations;
fp>>base_V;
fp>>base_VA;
base_I=base_VA/base_V;
fp>>no_DG;
fp>>no_tie;
getline(fp,comment);
getline(fp,comment);
Z_base=(base_V*base_V)/base_VA;
for(i=1;i<=no_lines;i++)
{
fp>>temp;
tran[i].isON=1;
fp>>tran[i].from;
fp>>tran[i].to;
fp>>tran[i].r;
fp>>tran[i].x;
z=complex<double>(tran[i].r,tran[i].x);
line_z[i]=z;
}
getline(fp,comment);
getline(fp,comment);

for(i=1;i<=no_buses;i++)
{
fp>>temp;
fp>>r;
fp>>x;
fp>>temp;
if(temp==1)
{
LD[i].alpha=1.51;
LD[i].beta=3.4;
}
if(temp==2)
{
LD[i].alpha=0.18;
LD[i].beta=6.0;
}
if(temp==3)
{
LD[i].alpha=0.92;
LD[i].beta=4.04;
}

LD[i].Pl0=r*1000/base_VA;
LD[i].Ql0=x*1000/base_VA;
}

getline(fp,comment);
getline(fp,comment);
for(i=1;i<=no_DG;i++)
{
fp>>temp;
fp>>r;
DG[temp].mp=r;
fp>>x;
DG[temp].nq=x;
fp>>x;
DG[temp].Qmax=x;
DG[temp].isDG=1;
fp>>r;
DG[temp].V0=r;
fp>>DG[temp].type;
}
getline(fp,comment);
getline(fp,comment);

for(i=no_lines+1;i<=(no_lines+no_tie);i++)
{
fp>>tran[i].from;
fp>>tran[i].to;
fp>>tran[i].r;
fp>>tran[i].x;
}
}

else
{cout<<"Unable to open file\n";
return 0;
}
fp.close();
tstep=1.0/1.0;
ntstep=1.0/tstep;
fp.open("loaddatai.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
getline(fp,comment);
getline(fp,comment);
l=1;
for(i=1;i<=no_hrs;i++)
{
fp>>r;
for(j=1;j<=ntstep;j++)
{
loaddatai[l]=r/100.0;
l++;
}
}
loaddatai[int(no_hrs*ntstep+1)]=loaddatai[1];
}
else
cout<<"Unable to open file\n";

fp.close();

fp.open("loaddatar.txt",ios::in);
if(fp.is_open())
{
getline(fp,comment);
getline(fp,comment);
getline(fp,comment);
l=1;
for(i=1;i<=no_hrs;i++)
{
fp>>r;
for(j=1;j<=ntstep;j++)
{
loaddatar[l]=r/100.0;
l++;
}
}
loaddatar[int(no_hrs*ntstep+1)]=loaddatar[1];
}
else
cout<<"Unable to open file\n";

fp.close();
fp.open("nrlfop.txt",ios::out);
fp<<std::fixed<<std::setprecision(4);
//									Change battery location here
location=8;
LBN=1;
Pdemand=0;
for(i=1;i<=no_buses;i++)
Pdemand+=LD[i].Pl0;

DG[33].a=4.356*1.17;
DG[33].b=0.345*1.17*1.5;
DG[33].c=0.0012*1.17/30.0;
DG[33].Pmax=(1500.0*1000.0)/base_VA;
DG[33].Smax=(1800.0*1000.0)/base_VA;

DG[18].a=4.526;
DG[18].b=0.3632;
DG[18].c=0.00017;
DG[18].Pmax=(1000.0*1000.0)/base_VA;
DG[18].Smax=(1000.0*1000.0)/base_VA;

DG[6].a=0.3312;
DG[6].b=0.0156;
DG[6].c=0.0002484;
DG[6].Pmax=(1800.0*1000.0)/base_VA;
DG[6].Smax=(1800.0*1000.0)/base_VA;

DG[LBN].a=0.4969*(2500.0/1790.0);
DG[LBN].b=0.0116;
DG[LBN].c=0.0001987*(1790.0/2500.0);
DG[LBN].Pmax=(2500.0*1000.0)/base_VA;

ncharge=sqrt(0.9);				//Charging effeiciency, %
ndis=ncharge;					//Discharging effeiciency, %
del=0.005;
delbat=100;
Pmax=0.2042;					//Power limit of the battery, kW
Emax=2.042;					//Storage capacity, kWh
Cbatt=0.3248;					//Cost of battery unit per day, $
socst=i2soc(soc2i(0.7));			//SoC at the start, %
cout<<"socst="<<socst<<"\n";
ndel=int(1.0/del)+1;
kmin=(Pmax*tstep)/(Emax*del*ndis);
kmax=(Pmax*tstep*ncharge)/(Emax*del);
cout<<"kmin="<<kmin<<"\tkmax="<<kmax<<"\n";
smin=0.4;					//Minimum limit of SoC, %
smax=1.0;					//Maximum limit of SoC, %
fp<<"socst="<<socst<<"\tM="<<(1.0/del)<<"\tlocation="<<location<<"\n";
socmin[1]=socst;
socmin1[1]=socst;
socmin2[1]=socst;
socmax[1]=socst;
socmax1[1]=socst;
socmax2[1]=socst;

socmin[int(no_hrs*ntstep+1)]=socst;
socmin1[int(no_hrs*ntstep+1)]=socst;
socmin2[int(no_hrs*ntstep+1)]=socst;
socmax[int(no_hrs*ntstep+1)]=socst;
socmax1[int(no_hrs*ntstep+1)]=socst;
socmax2[int(no_hrs*ntstep+1)]=socst;
for(i=0;i<=(1.0/del +50.0);i++)
for(j=0;j<=(no_hrs*ntstep+3);j++)
cost[i][j]=-999999999;
cost[soc2i(socst)][1]=0;
for(i=1;i<=(no_hrs*ntstep);i++)
{
socmin1[i+1]=max((socmin1[i]-kmin*del),smin);
socmin2[int(no_hrs*ntstep-i+1)]=max((socmin2[int(no_hrs*ntstep-i+2)]-kmax*del),smin);

socmax1[i+1]=min((socmax1[i]+kmax*del),smax);
socmax2[int(no_hrs*ntstep-i+1)]=min((socmax2[int(no_hrs*ntstep-i+2)]+kmin*del),smax);
}

for(i=1;i<=(no_hrs*ntstep);i++)
{
socmin[i+1]=max(socmin1[i+1],socmin2[i+1]);
socmax[i+1]=min(socmax1[i+1],socmax2[i+1]);
}
for(i=1;i<=(no_hrs*ntstep+1);i++)
cout<<socmin[i]<<"\t"<<socmax[i]<<"\n";
tcostold=9999999999999;
if(socst<smin)
{
cout<<"Incorrect data\n";
return 0;
}
n1=soc2i(0.0);
k=1;
Cq=0.01;
DG[LBN].isDG=1;
DG[6].isDG=1;
DG[33].isDG=0;
DG[18].isDG=1;
//									Comment following line to remove flow constraint
tran[6].limit=50;
LBN=1;
DG[LBN].type=1;
DG[6].type=6;
if(DG[LBN].isDG==0)
{
cout<<"LBN not a generator\n";
}
penalty[LBN]=1.0;
penalty[LBN+no_buses]=0;
DG[LBN].type=1;

for(i=1;i<=no_lines;i++)
{
bibc[i][tran[i].to-1]=1;

for(count=(tran[i].to-1);count<=no_buses;count++)
{
if(bibc[i][count])
{
no_to=findn(count+1);
for(j=1;j<=no_to;j++)
bibc[i][found[j]-1]=1;
}
}
}

for(k=1;k<=5000;k++)
{
if(nbatt==1700)
delbat=1;
cout<<nbatt<<"\n";
cost1=0;
m=10;
if(k==1)
{
fp<<"Solution of AC OPF without BESS\n";
fp<<"hour\t"<<"lambda"<<"\t"<<"DG[1].Pg"<<"\t"<<"DG[1].Qg"<<"\t"<<"DG[1].Sg"<<"\tpenalty[1]\t"<<"DG[6].Pg"<<"\t"<<"DG[6].Qg"<<"\t"<<"DG[6].Sg"<<"\t"<<"penalty[6]"<<"\t"<<"DG[33].Pg"<<"\t"<<"DG[33].Qg"<<"\t"<<"DG[33].Sg"<<"\t"<<"penalty[33]"<<"\n";
}
for(m=1;m<=no_hrs*ntstep;m++)
{
cout<<k<<"\t"<<m<<"\n";
n=2;
Pd=0;
for(i=1;i<=no_buses;i++)
if(i<23)
Pd+=LD[i].Pl0*loaddatar[m];
else
Pd+=LD[i].Pl0*loaddatai[m];
Pd+=battpower[m];
lambda=0.05;
for(i=1;i<=no_buses;i++)
{
muV[i]=0;
muIb[i]=0;
penalty[i]=1.0;
}
for(i=no_buses+2;i<=2*no_buses;i++)
penalty[i]=-0.0007;
e=0;
for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
{
DG[i].Qg=0;
}
//					First solve OPF without neglecting all constraints. Then include constraints one by one. 
do
{
tempsum1=0;
tempsum2=0;
for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
{
DG[i].Pg=(lambda - DG[i].b)/(2.0*DG[i].c);
tempsum2+=(1.0/(2.0*DG[i].c));
tempsum1+=DG[i].Pg;
}

lambda+=((Pd*base_VA/1000.0) - (tempsum1))/(tempsum2);
}while(abs((Pd*base_VA/1000.0) - (tempsum1))>(Pd*base_VA/1000000.0));
count=1;
do
{
tempsum1=0;
tempsum2=0;
qexit=0;

for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
{
double muVsump=0,muVsumq=0,mubsump=0,mubsumq=0;

DG[i].Pg=(lambda/penalty[i] -DG[i].Qper*Qpenalty[i]+muVsump -mubsump- DG[i].b)/(2.0*DG[i].c);

if(DG[i].Pg<(DG[i].Pmin*base_VA/1000.0))
{
DG[i].Pg=DG[i].Pmin*base_VA/1000.0;
}
if(DG[i].Pg>(DG[i].Pmax*base_VA/1000.0))
{
DG[i].Pg=DG[i].Pmax*base_VA/1000.0;
}

tempsum1+=DG[i].Pg;
tempsum2+=DG[i].Pmin;
DG[i].Pg*=(1000.0/base_VA);

DG[i].Qg+=9000.0*(1000.0/base_VA)*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses])+muVsumq-mubsumq);
if(DG[i].Qg>sqrt(DG[i].Smax*DG[i].Smax-DG[i].Pg*DG[i].Pg))
{
DG[i].Qg=sqrt(DG[i].Smax*DG[i].Smax-DG[i].Pg*DG[i].Pg);
}
else
if((qexit<(9000.0*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses]))))&&(i!=1))
qexit=(9000.0*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses])+muVsumq-mubsumq));
}

if(tempsum2>Pd)
{
cout<<"Sorry not feasible\n";
break;
}

islandnrlf(m,battpower[m]);
for(i=1;i<=no_lines;i++)
{
tempsum=complex<double>(0,0);
for(j=2;j<=no_buses;j++)
tempsum+=I[j]*complex<double>(bibc[i][j-1],0);
tran[i].flow=tempsum;
}
lambda+=0.0001*(-tempsum1+Ploss+Pd*base_VA/1000.0);
for(i=1;i<=(2*no_buses);i++)
Jlbn[i]=J1[2*no_buses-1][i];
for(i=1;i<=(2*no_buses);i++)
Qlbn[i]=J1[2*no_buses][i];

invshipley(J1,(2*no_buses-2));

for(i=1;i<=(2*no_buses-2);i++)
{
temp1=0;
for(j=1;j<=(2*no_buses-2);j++)
temp1+=Jlbn[j]*J1[j][i];
if((i+1)<=no_buses)
penalty[i+1]=(-1.0/temp1);
else
penalty[i+2]=temp1;
temp1=0;
for(j=1;j<=(2*no_buses-2);j++)
temp1+=Qlbn[j]*J1[j][i];
if((i+1)<=no_buses)
Qpenalty[i+1]=temp1;
else
Qpenalty[i+2]=temp1;
}
}while(!(abs((Pd*base_VA/1000.0) + Ploss - (tempsum1))<(1.0)&&(qexit<5.0)));
do
{
tempsum1=0;
tempsum2=0;
qexit=0;

for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
{
double muVsump=0,muVsumq=0,mubsump=0,mubsumq=0;

for(int count=1;count<=no_buses;count++)
{
if(abs(V[count])<Vmin)
if(count==1)
{
muVsump+=muV[count]*J[2*no_buses][2*no_buses-1]*abs(V[count]);
muVsumq+=muV[count]*J[2*no_buses][2*no_buses]*abs(V[count]);
}
else
{
muVsump+=muV[count]*J[count+no_buses-2][count-1]*abs(V[count]);
muVsumq+=muV[count]*J[count+no_buses-2][count+no_buses-2]*abs(V[count]);
}

}
for(count=1;count<=no_lines;count++)
if(muIb[count])
{
int iter;
double sumr=0,sumi=0;
for(iter=2;iter<=no_lines;iter++)
{
sumr+=bibc[count][iter-1]*(-imag(I[iter])*J[iter-1][i-1]-real(I[iter])*J[no_buses-2+iter][i-1]);
sumr+=bibc[count][i-1]*(cos(arg(V[i]))/abs(V[i]));
sumi+=bibc[count][iter-1]*(real(I[iter])*J[iter-1][i-1]-imag(I[iter])*J[no_buses-2+iter][i-1]);
sumi+=bibc[count][i-1]*(sin(arg(V[i]))/abs(V[i]));
}
mubsump+=muIb[count]*(1.0/abs(tran[count].flow))*(real(tran[count].flow)*sumr +imag(tran[count].flow)*sumi);
sumr=0;
sumi=0;
for(iter=2;iter<=no_lines;iter++)
{
sumr+=bibc[count][iter-1]*(-imag(I[iter])*J[iter-1][i+no_buses-2]-real(I[iter])*J[no_buses-2+iter][i+no_buses-2]);
sumr+=bibc[count][i-1]*(sin(arg(V[i]))/abs(V[i]));
sumi+=bibc[count][iter-1]*(real(I[iter])*J[iter-1][i+no_buses-2]-imag(I[iter])*J[no_buses-2+iter][i+no_buses-2]);
sumi-=bibc[count][i-1]*(cos(arg(V[i]))/abs(V[i]));
}
mubsumq+=muIb[count]*(1.0/abs(tran[count].flow))*(real(tran[count].flow)*sumr +imag(tran[count].flow)*sumi);
}
DG[i].Pg=(lambda/penalty[i] -DG[i].Qper*Qpenalty[i]+muVsump -mubsump- DG[i].b)/(2.0*DG[i].c);
if(DG[i].Pg<(DG[i].Pmin*base_VA/1000.0))
{
DG[i].Pg=DG[i].Pmin*base_VA/1000.0;
}
if(DG[i].Pg>(DG[i].Pmax*base_VA/1000.0))
{
DG[i].Pg=DG[i].Pmax*base_VA/1000.0;
}

tempsum1+=DG[i].Pg;
tempsum2+=DG[i].Pmin;
DG[i].Pg*=(1000.0/base_VA);
DG[i].Qg+=1000.0*(1000.0/base_VA)*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses])+muVsumq-mubsumq);
if(DG[i].Qg>sqrt(DG[i].Smax*DG[i].Smax-DG[i].Pg*DG[i].Pg))
{
DG[i].Qg=sqrt(DG[i].Smax*DG[i].Smax-DG[i].Pg*DG[i].Pg);
}
else
if((qexit<(1000.0*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses]))))&&(i!=1))
qexit=(1000.0*(-lambda*penalty[i+no_buses]-DG[i].Qper*(1.0+Qpenalty[i+no_buses])+muVsumq-mubsumq));
}
if(tempsum2>Pd)
{
cout<<"Sorry not feasible\n";
break;
}
islandnrlf(m,battpower[m]);
for(i=1;i<=no_lines;i++)
{
tempsum=complex<double>(0,0);
for(j=2;j<=no_buses;j++)
tempsum+=I[j]*complex<double>(bibc[i][j-1],0);
tran[i].flow=tempsum;
}
if((m==10))
lstep=0.00005;
else
lstep=0.0001;
lambda+=lstep*(-tempsum1+Ploss+Pd*base_VA/1000.0);
for(i=1;i<=no_buses;i++)
if(Vmin>abs(V[i]))
{
muV[i]+=0.01*(Vmin-abs(V[i]));
}
for(i=1;i<=no_lines;i++)
if(abs(tran[i].flow)*base_I>tran[i].limit)
{
muIb[i]+=0.00005*(abs(tran[i].flow)*base_I-tran[i].limit);
}
for(i=1;i<=(2*no_buses);i++)
Jlbn[i]=J1[2*no_buses-1][i];
for(i=1;i<=(2*no_buses);i++)
Qlbn[i]=J1[2*no_buses][i];
invshipley(J1,(2*no_buses-2));

for(i=1;i<=(2*no_buses-2);i++)
{
temp1=0;
for(j=1;j<=(2*no_buses-2);j++)
temp1+=Jlbn[j]*J1[j][i];
if((i+1)<=no_buses)
penalty[i+1]=(-1.0/temp1);
else
penalty[i+2]=temp1;
temp1=0;
for(j=1;j<=(2*no_buses-2);j++)
temp1+=Qlbn[j]*J1[j][i];
if((i+1)<=no_buses)
Qpenalty[i+1]=temp1;
else
Qpenalty[i+2]=temp1;
}
}while(!(abs((Pd*base_VA/1000.0) + Ploss - (tempsum1))<(1.0)&&(qexit<5.0)&&(Vmin-min(V,no_buses)<0.005)&&(abs(tran[6].flow)*base_I-tran[6].limit<0.1)));

LMPE[1][m]=lambda;
LMP[1][m]=lambda;
for(i=2;i<=no_buses;i++)
{
double muVsump=0,muVsumq=0,mubsump=0,mubsumq=0;

for(int count=1;count<=no_buses;count++)
{
if(count==1)
{
muVsump+=muV[count]*J[2*no_buses][2*no_buses-1]*abs(V[count]);
muVsumq+=muV[count]*J[2*no_buses][2*no_buses]*abs(V[count]);
}
else
{
muVsump+=muV[count]*J[count+no_buses-2][count-1]*abs(V[count]);
muVsumq+=muV[count]*J[count+no_buses-2][count+no_buses-2]*abs(V[count]);
}

}
for(count=1;count<=no_lines;count++)
if(muIb[count])
{
int iter;
double sumr=0,sumi=0;
for(iter=2;iter<=no_lines;iter++)
{
sumr+=bibc[count][iter-1]*(-imag(I[iter])*J[iter-1][i-1]-real(I[iter])*J[no_buses-2+iter][i-1]);
sumr+=bibc[count][i-1]*(cos(arg(V[i]))/abs(V[i]));
sumi+=bibc[count][iter-1]*(real(I[iter])*J[iter-1][i-1]-imag(I[iter])*J[no_buses-2+iter][i-1]);
sumi+=bibc[count][i-1]*(sin(arg(V[i]))/abs(V[i]));
}
mubsump+=muIb[count]*(1.0/abs(tran[count].flow))*(real(tran[count].flow)*sumr +imag(tran[count].flow)*sumi);
sumr=0;
sumi=0;
for(iter=2;iter<=no_lines;iter++)
{
sumr+=bibc[count][iter-1]*(-imag(I[iter])*J[iter-1][i+no_buses-2]-real(I[iter])*J[no_buses-2+iter][i+no_buses-2]);
sumr+=bibc[count][i-1]*(sin(arg(V[i]))/abs(V[i]));
sumi+=bibc[count][iter-1]*(real(I[iter])*J[iter-1][i+no_buses-2]-imag(I[iter])*J[no_buses-2+iter][i+no_buses-2]);
sumi-=bibc[count][i-1]*(cos(arg(V[i]))/abs(V[i]));
}
mubsumq+=muIb[count]*(1.0/abs(tran[count].flow))*(real(tran[count].flow)*sumr +imag(tran[count].flow)*sumi);
}
LMPV[i][m]=muVsump;
LMPI[i][m]=-mubsump;
LMPloss[i][m]=lambda*(1.0-penalty[i])/penalty[i];
LMP[i][m]=lambda+LMPV[i][m]+LMPI[i][m]+LMPloss[i][m];
LMPE[i][m]=lambda;
}
islandnrlf(m,battpower[m]);
sortdescend(V,no_buses);
for(i=1;i<=no_buses;i++)
if(DG[int(Vsort[i][2])].isDG==0)
{
if(Vsort[i][1]>=Vmax)
cout<<m<<"\tVmax=V"<<Vsort[i][2]<<"="<<Vsort[i][1]<<"\n";
break;
}
for(i=no_buses;i>=1;i--)
{
if(Vsort[i][1]<=Vmin)
cout<<m<<"\tVmin=V"<<Vsort[i][2]<<"="<<Vsort[i][1]<<"\t"<<muV[int(Vsort[i][2])]<<"\n";
break;
}
for(i=1;i<=no_lines;i++)
{
tempsum=complex<double>(0,0);
for(j=2;j<=no_buses;j++)
tempsum+=I[j]*complex<double>(bibc[i][j-1],0);
tran[i].flow=tempsum;
if(abs(tran[i].flow)*base_I>=tran[i].limit)
cout<<m<<"\t"<<i<<"\t"<<abs(tran[i].flow)*base_I<<"\n";
}
if(k==1)
{
fp<<m<<"\t"<<lambda<<"\t";
for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
fp<<DG[i].Pg*base_VA/1000.0<<"\t"<<DG[i].Qg*base_VA/1000.0<<"\t"<<sqrt(DG[i].Qg*DG[i].Qg+DG[i].Pg*DG[i].Pg)*base_VA/1000.0<<"\t"<<penalty[i]<<"\t";
fp<<abs(tran[6].flow)*base_I<<"\n";
}
price[m]=LMP[location][m];
for(i=1;i<=no_buses;i++)
if(DG[i].isDG==1)
{
cost1+=(DG[i].a+DG[i].b*DG[i].Pg*base_VA/1000.0 + DG[i].c*(DG[i].Pg*base_VA/1000.0)*(DG[i].Pg*base_VA/1000.0))*tstep;
cost1+=DG[i].Qper*abs(DG[i].Qg)*base_VA/1000.0;
}
}
cout<<"Total cost="<<(cost1 +(k-1)*Cbatt)<<"$/day\n";
if(k==1)
{
fp<<"Total DG cost="<<cost1<<"$/day\n";
fp<<"Total battery cost="<<((k-1)*Cbatt)<<"$/day\n";
fp<<"Total cost="<<(cost1 +(k-1)*Cbatt)<<"$/day\n";

fp<<"LMPtotal (Node no, hour, value)\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<"("<<i<<","<<j<<","<<LMP[i][j]<<")"<<"\t";
fp<<"\n";
}
fp<<"\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMP[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPEnergy\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPE[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPloss\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPloss[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPVoltage\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPV[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPflow\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPI[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
}
for(j=2;j<=(no_hrs*ntstep+1);j++)
for(i=soc2i(1.0);i<=n1;i++)
{
if((i2soc(i)-socmax[j]>0.000001)||(socmin[j]-i2soc(i)>0.000001))
{
cost[i][j]=-999999999;
}
else
{
l=1;
for(count=-kmin;count<=kmax;count++)
{
if((i+count)<0)
{
a[l]=-999999999;
l++;
continue;
}
if(count<0)
a[l]=cost[i+count][j-1]-del*count*Emax*ndis*price[j-1];
else
a[l]=cost[i+count][j-1]-del*count*Emax*price[j-1]/ncharge;
l++;
}
cost[i][j]=max(a,(kmin+kmax+1));
}}
i=soc2i(socst);
temp1=0;
soc=socst;
charge[int(no_hrs*ntstep+1)]=0;
for(j=(no_hrs*ntstep+1);j>=1;j--)
{
l=1;
for(count=-kmin;count<=kmax;count++)
{
if((i+count)<0)
{
a[l]=-999999999;
l++;
continue;
}
if(count<0)
{
if(cost[i][j]==(cost[i+count][j-1] -del*count*Emax*ndis*price[j-1]))
{
charge[j]=count;
break;
}
}
else
if(cost[i][j]==(cost[i+count][j-1] -del*count*Emax*price[j-1]/ncharge))
{
charge[j]=count;
break;
}
l++;
}
i+=charge[j];
if(charge[j]>0)
temp1-=charge[j]*price[j-1]*del*Emax/ncharge;
else
temp1-=charge[j]*price[j-1]*del*Emax*ndis;
soc-=charge[j]*del;
}

if((cost1 +nbatt*Cbatt)<tcostold)
{
tcostold=(cost1 +nbatt*Cbatt);
nbatt+=delbat;
for(m=1;m<=(no_hrs*ntstep);m++)
{
if(charge[m+1]>0)
battpower[m]+=(delbat*charge[m+1]*del*Emax*1000.0/(base_VA*ncharge))/tstep;
else
battpower[m]+=(ndis*delbat*charge[m+1]*del*Emax*1000.0/(base_VA))/tstep;
}
}
else
{
cout<<"stopped after "<<k<<"\n";
break;
}
}
socT[1]=socst;
for(m=1;m<=(no_hrs*ntstep);m++)
{
if(battpower[m]>0)
socT[m+1]=socT[m]+ncharge*battpower[m]*tstep*base_VA/(1000.0*nbatt*Emax);
else
{
socT[m+1]=socT[m]+battpower[m]*tstep*base_VA/(1000.0*nbatt*Emax*ndis);
}
}
fp<<"\n\nOptimal # of batteries="<<nbatt<<"\n"<<"n="<<(ncharge*ndis)<<"\tsocmin="<<smin<<"\tCbatt="<<Cbatt<<"\n";
fp<<"Total DG cost="<<cost1<<"$/day\n";
fp<<"Total battery cost="<<(nbatt*Cbatt)<<"$/day\n";
fp<<"Total cost="<<(cost1 +nbatt*Cbatt)<<"$/day\n";
fp<<"\nOptimal battery power (Pmax="<<Pmax<<" kW Emax="<<Emax<<" kWh)\n";
fp<<"hour\tpower(kW)\tpower per battery\t\tin pu\t\tSOC\n";
for(m=1;m<=(no_hrs*ntstep+1);m++)
{
fp<<m<<"\t"<<battpower[m]*base_VA/1000.0<<"\t\t"<<battpower[m]*base_VA/(1000.0*nbatt)<<"\t\t\t"<<battpower[m]*base_VA/(1000.0*nbatt*Pmax)<<"\t\t"<<socT[m]<<"\n";
}
fp<<"LMP with BESS\n";
fp<<"LMPtotal\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<"("<<i<<","<<j<<","<<LMP[i][j]<<")"<<"\t";
fp<<"\n";
}
fp<<"\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMP[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPEnergy\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPE[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPloss\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPloss[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPVoltage\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPV[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp<<"LMPflow\n";
for(i=1;i<=no_buses;i++)
{
for(j=1;j<=no_hrs;j++)
fp<<LMPI[i][j]<<"\t";
fp<<"\n";
}
fp<<"\n";
fp.close();
return 0;
}
