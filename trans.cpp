#include "alp.h"
// transmat fuancion ma(neV),e(KeV),ne(cm^-3) gag(10^-11 GeV^-1),b(10^-6 G) deta (Kpc)-1

cmat transmat(double ma,double ga,double e,double b,double fi,double ne,double l)
{
    double dag=1.52e-2*(ga)*(b/1.0e-3)*1.0e-3;
    double da=-7.8e-3*ma*ma/1.0e-8/(e/100)*1.0e-3;
    double dpl=-1.1e-4/(e/100)*ne/1.0e-7*1.0e-3;
    double dqed=4.1e-16*(e/100)*(b/1.0e-3)*(b/1.0e-3)*1.0e-3;
    double dgg=8.0e-2*(e/1.0e9)*1.0e-3;
    //double dt=dpl,dp=dpl;
    //double dt=dpl+2.0*dqed+dgg, dp=dpl+3.5*dqed+dgg;
    double dt=0, dp=0; // for comparing with Libanov2020
    double theta=0.5*atan(2*dag/(dp-da));
    if(theta<0.0) theta=theta+pi/2;
    double d1=dt;
    double d2=(dp+da)/2+0.5*sqrt((da-dp)*(da-dp)+4*dag*dag);
    double d3=(dp+da)/2-0.5*sqrt((da-dp)*(da-dp)+4*dag*dag);

    double sfi=sin(fi),cfi=cos(fi);
    double sfi2=sin(2*fi);
    double sth=sin(theta),cth=cos(theta);
    double sth2=sin(2*theta);
    mat ta(0.0);
    ta.m[0][0]=pow(cfi,2);
    ta.m[0][1]=-sfi2/2;
    ta.m[1][0]=-sfi2/2;
    ta.m[1][1]=pow(sfi,2);
    mat tb(0.0);
    tb.m[0][0]=pow(sfi*cth,2);
    tb.m[0][1]=sfi2/2*pow(cth,2);
    tb.m[0][2]=sth2/2*sfi;
    tb.m[1][0]=tb.m[0][1];
    tb.m[1][1]=pow(cfi*cth,2);
    tb.m[1][2]=sth2/2*cfi;
    tb.m[2][0]=tb.m[0][2];
    tb.m[2][1]=tb.m[1][2];
    tb.m[2][2]=pow(sth,2);
    mat tc(0.0);
    tc.m[0][0]=pow(sfi*sth,2);
    tc.m[0][1]=sfi2/2*pow(sth,2);
    tc.m[0][2]=-sth2/2*sfi;
    tc.m[1][0]=tc.m[0][1];
    tc.m[1][1]=pow(sth*cfi,2);
    tc.m[1][2]=-sth2/2*cfi;
    tc.m[2][0]=tc.m[0][2];
    tc.m[2][1]=tc.m[1][2];
    tc.m[2][2]=pow(cth,2);
    //ta.print();
    //tb.print();
    //tc.print();
    //std::cout<<"d1:"<<d1<<"  "<<d2<<"  "<<d3<<std::endl;
    //if(b>6.3e10) std::cout<<"d1:"<<b<<"  "<<da<<"   "<<"  "<<dpl<<"  "<<std::endl;
    cmat mattmp=cmat(ta)*cexp(cplex(0,-l*d1))+cmat(tb)*cexp(cplex(0,-l*d2))+cmat(tc)*cexp(cplex(0,-l*d3));
 //cmat mattmp=cmat(mat(0.0));
    //mattmp.print();
	return mattmp;
}

