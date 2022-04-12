/*Problema de biofotólise indireta - Tales*/

#include<stdio.h>
#include<math.h>
#include<time.h>

#define ikmax 101*powf(10,-6)
#define ik 363*powf(10,-6)
#define Io 200.0
#define A1 3.2*powf(10,5)
#define A2 1.44*powf(10,12)
#define Ea1 3.75*powf(10,4)
#define Ea2 7.85*powf(10,4)
#define absor 7*powf(10,-3)
#define mumax1 0.1
#define ks1 49.9
#define yx1r 1.11
#define tmax 180.0
#define h 0.01
#define L 1.0
#define mumax2 0.210
#define yxp 0.1333
#define mumax3 7.65417*powf(10,-6)
#define yxh2 0.10
#define taux 2
#define thidro 41


double Mu1(double R)
{
    double m1=mumax1*R/(R+ks1);
    return(m1);
}
 double Mu2()
{
    double phi,m2,Im;
    phi = ikmax*Io/(ik+Io);
    Im=Io*(1-exp(absor*L))/absor*L;
    m2=mumax2*Im/(phi+Im);
    return(m2);
}
double Mu3(double S)
{
    double m3=mumax3*S/(S+ks1);
    return(m3);
}
double *F1(double R,double S,double X,double *edo)
{
    double a, mu;
    mu=Mu1(R);
    a=h*mu*X;
    edo[0]=-a*yx1r;
    edo[1]=-edo[0];
    edo[2]=a;
}
double *F2(double S2,double X2,double *edo2)
{
    double b, mu2;
    mu2=Mu2();
    b=h*mu2*X2;
    edo2[0]=b*yxh2;
    edo2[1]=b;

}
double *F3(double P,double S,double X2,double *edo3)
{
    double c, mu3;
    mu3=Mu3(S);
    c=h*mu3*X2;
    edo3[0]=-c*yxp;
    edo3[1]=12*c;
}
double RK(double R,double S,double X, double t)
{
    double tzero, e=1.0/6.0,k[4],l[4],m[4],U1[3],U2[3],U3[3],U4[3];
    tzero=t;
    FILE *rk=fopen("RungeKutta.txt","w");
    while(tzero<=taux)
    {
        fprintf(rk,"%f %f %f %f\n",tzero,R,S,X);
        F1(R,S,X,U1);
        k[0]=U1[0];
        l[0]=U1[1];
        m[0]=U1[2];
        F1(R+0.5*k[0],S+0.5*l[0],X+0.5*m[0],U2);
        k[1]=U2[0];
        l[1]=U2[1];
        m[1]=U2[2];
        F1(R+0.5*k[1],S+0.5*l[1],X+0.5*m[1],U3);
        k[2]=U3[0];
        l[2]=U3[1];
        m[2]=U3[2];
        F1(R+k[2],S+l[2],X+m[2],U4);
        k[3]=U4[0];
        l[3]=U4[1];
        m[3]=U4[2];
        /*Atualização*/
        R=R+e*(k[0]+2*(k[1]+k[2])+k[3]);
        S=S+e*(l[0]+2*(l[1]+l[2])+l[3]);
        X=X+e*(m[0]+2*(m[1]+m[2])+m[3]);
        tzero=tzero+h;
    }
    fclose(rk);
}
double RK2(double S2,double X2, double t)
{

    double tzero,e=1.0/6.0,k2[4],l2[4],m2[4],U12[2],U22[2],U32[2],U42[2];
    tzero=t;
    FILE *rk2=fopen("RungeKutta2.txt","w");
    while(tzero<=taux)
    {
        fprintf(rk2,"%f %f %f\n ",tzero,S2,X2);
        F2(S2,X2,U12);
        k2[0]=U12[0];
        l2[0]=U12[1];
        F2(S2+0.5*k2[0],X2+0.5*l2[0],U22);
        k2[1]=U22[0];
        l2[1]=U22[1];
        F2(S2+0.5*k2[1],X2+0.5*l2[1],U32);
        k2[2]=U32[0];
        l2[2]=U32[1];
        F2(S2+k2[2],X2+l2[2],U42);
        k2[3]=U42[0];
        l2[3]=U42[1];
        /*Atualização*/
        X2=X2+e*(l2[0]+2*(l2[1]+l2[2])+l2[3]);
        S2=S2+e*(k2[0]+2*(k2[1]+k2[2])+k2[3]);
        tzero=tzero+h;
    }
    fclose(rk2);
}

 double RK3(double P,double S3,double X2,double t)
{
    double tzero,e=1.0/6.0,k3[4],l3[4],m3[4],U13[3],U23[3],U33[3],U43[3];
    tzero=t;
    FILE *rk3=fopen("RungeKutta3.txt","w");
    while(tzero<=taux)
    {
        fprintf(rk3,"%f %f %f %f\n",tzero,X2,P,S3);
        F3(P,S3,X2,U13);
        k3[0]=U13[0];
        l3[0]=U13[1];
        m3[0]=U13[2];
        F3(P+0.5*k3[0],S3+0.5*l3[0],X2+0.5*m3[0],U23);
        k3[1]=U23[0];
        l3[1]=U23[1];
        m3[1]=U23[2];
        F3(P+0.5*k3[1],S3+0.5*l3[1],X2+0.5*m3[1],U33);
        k3[2]=U33[0];
        l3[2]=U33[1];
        m3[2]=U33[2];
        F3(P+k3[2],S3+l3[2],X2+m3[2],U43);
        k3[3]=U43[0];
        l3[3]=U43[1];
        m3[3]=U43[2];
        /*Atualização*/
        P=P+e*(k3[0]+2*(k3[1]+k3[2])+k3[3]);
        S3=S3+e*(l3[0]+2*(l3[1]+l3[2])+l3[3]);
        X2=X2+e*(m3[0]+2*(m3[1]+m3[2])+m3[3]);
        tzero=tzero+h;
    }
    fclose(rk3);
}

main()
{
    int i,j,n,taux1,N2,N,T;
    double R,S,X,P,v,mu,m,m1,m2,Mu2,phi,Im,mu2,S2,X2,Mu3,m3,M[N][4],F2[N2][3],F3[N2][4],SF0;
    RK(5,0,0.5,0);
    X2=0.5;
    v=1/h;
    N=taux/v;
    P=0;
    N2=12/v;

    FILE *rk=fopen("RungeKutta.txt","r");

    for(i=0;i<N;i++)
    {
        for(j=0;j<4;j++)
        {
            fscanf(rk,"%lf",&M[i][j]);
        }
        fscanf(rk,"\n");
    }
    fclose(rk);
    SF0=M[N][1];
    printf("%f\n",SF0);
    taux1=0;
    T=0;
    while(T<360)
    {
        printf("%d\n",taux1);
        if(taux1==0)
        {
            RK2(SF0,X2,taux1);
            FILE *rk2=fopen("RungeKutta2.txt","r");
            for(i=0;i<N2;i++)
            {
                for(j=0;j<3;j++)
                {
                    fscanf(rk2,"%lf",&F2[i][j]);
                }
                fscanf(rk2,"\n");
            }
            fclose(rk2);
            SF0=F2[N2][1];
            X2=F2[N2][2];
            taux1=12;
        }
        else
        {
            RK3(SF0,P,X2,taux);
            FILE *rk3=fopen("RungeKutta3.txt","r");
            for(i=0;i<N2;i++)
            {
                for(j=0;j<3;j++)
                {
                    fscanf(rk3,"%lf",&F3[i][j]);
                    printf("%f ",M[i][j]);
                }
                fscanf(rk3,"\n");
                printf("\n");
            }
            fclose(rk3);
            SF0=F3[N2][3];
            X2=F3[N2][1];
            taux1=0;
        }
    T=T+12;
    }
}