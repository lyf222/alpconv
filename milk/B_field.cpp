#include"alp.h"

//Jansson&Farrar 2012 model
// 发现z严格=0时候 即在银盘上l=0时 会出现NaN 后面再解决这个问题 暂时先避免l=0.0

inline double lfun(double z,double h,double w)
{
    return 1.0/(1.0+exp(-2.0*(fabs(z)-h)/w));
}
vec cord(double ra,double dec)
{
    double ra0=192.859481/180*pi;
    double d0=27.12825118/180*pi;
    double l0=122.931918568/180*pi;
    double bb=asin(sin(dec)*sin(d0)+cos(dec)*cos(d0)*cos(ra-ra0));
    double y=cos(dec)*sin(ra-ra0);
    double x=sin(dec)*cos(d0)-cos(dec)*sin(d0)*cos(ra-ra0);
    double ll=atan(y/x);
    if(x>0)
    {
        if(y<0) ll=ll+2*pi;
    }
    else
    {
        ll=ll+pi;
    }
    ll=l0-ll;
    if(ll<0) ll=ll+2*pi;
    vec cordb_xy(0.0);
    std::cout<<ra*180/pi<<"  "<<dec*180/pi<<std::endl;
    std::cout<<ll*180/pi<<"  "<<bb*180/pi<<std::endl;
    cordb_xy.v[2]=cos(bb);
    cordb_xy.v[1]=sin(bb)*sin(ll);
    cordb_xy.v[0]=sin(bb)*cos(ll);
    //std::cout<<"ra="<<ra*180/pi<<"  dec"<<dec*180/pi<<std::endl;
    //std::cout<<"ll="<<ll*180/pi<<"  bb"<<bb*180/pi<<std::endl;
    return cordb_xy;
}

// coherent magnetic field in Milky Way and turbulent magnetic field in gamma ray source
vec b_fieldm(vec cordxy)
{
    double x=cordxy.v[0];
    double y=cordxy.v[1];
    double z=cordxy.v[2];
    double r=sqrt(x*x+y*y);
    if(r>20) return vec(0);
    if(cordxy.mag()<1) return vec(0);
    double fi=atan(y/x);
    if(x>0)
    {fi=fi+pi;}
    else if(y>0)
    {
        fi=fi+2*pi;
    }

    if(x==0)
    {
        if(y<0)
        {
            fi=pi/2;
        }
        else
        {
            fi=3*pi/2;
        }
    }
    vec bdisk(0.0),bhalo(0.0),bxx(0.0);
    // disk component
    if(r>3.0){
        if(r<5.0){
            bdisk.v[1]=bring;
        }
        else{
            double tmp[24];
            for(int i=0;i<3;i++)
            {
                double tmpfi=fi+(i-1)*2*pi;
                for(int j=0;j<8;j++)
                {
                    tmp[i*8+j]=sdisk[j]*exp(tmpfi*tan(openi));
                }
            }
            int id=0;
            for(int i=0;i<24;i++)
            {
                if(r<tmp[i])
                {
                    id=i;
                    break;
                }
            }
            //std::cout<<id%8<<" ";
            id=id%8;
            bdisk.v[0]=bd[id]/r*5*sin(openi);
            bdisk.v[1]=bd[id]/r*5*cos(openi);
        }}
    bdisk=bdisk*(1-lfun(z,hdisk,wdisk));

    //halo component
    double factor=exp(-fabs(z)/z0)*lfun(z,hdisk,wdisk);
    if(z>0)
    {
        bhalo.v[1]=factor*bn*(1-lfun(r,rn,wh));
    }
    if(z<0)
    {
        bhalo.v[1]=factor*bs*(1-lfun(r,rs,wh));
    }

    // X-field component
    double rp=r-fabs(z)/tan(theta0);
    double bth=0.0,theta=theta0;
    if(rp>=rxc)
    //if(cordxy.mag()>rxc)
    {
        bth=bx*exp(-rp/rx)*rp/r;
    }
    else
    {
        rp=r*rxc/(rxc+fabs(z)/tan(theta0));
        theta=atan(fabs(z)/(r-rp));
        bth=bx*exp(-rp/rx)*(rp/r)*(rp/r);
    }
    bxx.v[2]=bth*sin(theta);
    if(z>0) bxx.v[0]=bth*cos(theta);
    if(z<0) bxx.v[0]=-bth*cos(theta);
    vec vectmp=bdisk+bhalo+bxx;
    //vec vectmp=bhalo+bxx;
    //vec vectmp=bxx;
    //vec vectmp=bdisk;
    vec vecxy(0.0);
    vecxy.v[0]=x/r*vectmp.v[0]-y/r*vectmp.v[1];
    vecxy.v[1]=y/r*vectmp.v[0]+x/r*vectmp.v[1];
    vecxy.v[2]=vectmp.v[2];

    return vecxy;
}
double fq1(double q,double s1,double s2)
{
    double tmp=(pow(s2,q)-pow(s1,q))/q;
    if(abs(q) < 1.e-30) return log(s2/s1);
    return tmp;
}
double fq(double q,double k,double kl,double kh)
{
    if(k>=kh) return 0.0;
    double kk=k>kl?k:kl;
    if(k>kl)
    {
        kk=k;
    }
    else
    {
        kk=kl;
    }
    return (fq1(q+2,kk,kh)+k*k*fq1(q,kk,kh))/fq1(q+3,kl,kh);
}
field b_fields(field& random,field& random1)
{
    double a=sigamb*sigamb/8/pi*fq1(q+3,kl,kh);
    field b(nl,0.0);
    for(int i=0;i<nk;i++)
    {
        double k=lmin*pow(stepk,i);
        double deltk=k*(stepk-1);
        double et=pi*sigamb*sigamb/4*fq(q,k,kl,kh);
        double factor=sqrt(2*et*deltk/pi*log(1.0/random.u[i]));
        double phase=2*pi*random.v[i];
        double factor1=sqrt(2*et*deltk/pi*log(1.0/random1.u[i]));
        double phase1=2*pi*random1.v[i];

        for(int j=0;j<nl;j++)
        {
            double x3=stepl*j;
            b.u[j]+=factor*cos(k*x3+phase);
            b.v[j]+=factor1*cos(k*x3+phase1);
        }
    }
    double f0=edens(0.0);
	for(int j=0;j<nl;j++)
    {
        double x3=stepl*j;
        double f=pow(edens(x3)/f0,yita);
        b.u[j]*=f;
        b.v[j]*=f;
    }
    return b;
}
