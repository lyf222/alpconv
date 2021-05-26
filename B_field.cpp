#include"alp.h"
#include"utils.h"
#include <fstream>
using std::cout;
using std::endl;
using std::vector;
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
    //std::cout<<q<<"  "<<tmp<<std::endl;
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


/* contents copy from alp.h
//parameter for magnetic field in 星系团磁场
//const double kl=2*pi/35,kh=2*pi/0.7,q=-2.8,rmax=370,sigamb=10.0,yita=0.5;
const double kl=2*pi/35,kh=2*pi/0.7,q=-2.8,rmax=500.,sigamb=10.0,yita=0.5;

const int nk=int(20*(log10(kh/kl)+3));  // nk is the total number of spacings in k; where k is the wave number
const double stepk=pow(10,(log10(kh/kl)+3)/nk);

const int nl= (rmax*kh)/4; // l is the integration step length, thus nl is the number of step
const double lmin=kl*1.0e-3,stepl=rmax/nl;
*/

field b_fields(field& random,field& random1)
{
    //double a=sigamb*sigamb/8/pi*fq1(q+3,kl,kh);
    field b(nl,0.0);
    //std::cout<<nl<<"abcdefg"<<std::endl;
    int nkk = 10.*(log10(kh)-log10(lmin))*(log10(kh)-log10(lmin))+1; //Meyer程序中的取nk方式
    //std::cout<<nkk<<"abcdefg"<<nl<<"nl"<<rmax<<"rmax"<<std::endl;
    //10**(1+((4-1.)/6.)*6)

    for(int j=0;j<nl;j++)
    {
        double x3=stepl*(j+0.5); //stepl*j;
        //double x3=0.5; //for debug
        for(int i=0;i<nk-1;i++)
        {
            //double k=lmin*pow(stepk,i); //Z
            double stepkk=(log10(kh)-log10(lmin))/(nk-1);
            double kk=log10(lmin)+i*stepkk; //me
            double kk1=log10(lmin)+(i+1)*stepkk; //me
            double deltkk=pow(10,kk1)-pow(10,kk); //me
            //double deltk=k*(stepk-1); //Z
            double k=pow(10,kk);
            double et=pi*sigamb*sigamb/4*fq(q,k,kl,kh); //Z
            double factor=sqrt(2*et*deltkk/pi*log(1.0/random.u[i])); //Z
            double phase=2*pi*random.v[i]; //Z
            double factor1=sqrt(2*et*deltkk/pi*log(1.0/random.u[i])); //Z
            double phase1=2*pi*random.v[i]; //Z
            //double factor1=sqrt(2*et*deltkk/pi*log(1.0/random1.u[i])); //Z
            //double phase1=2*pi*random1.v[i]; //Z

            b.u[j]+=factor*cos(k*x3+phase); //Z
            b.v[j]+=factor1*cos(k*x3+phase1); //Z
            //std::cout<<'i'<<i<<"kk"<<pow(10,kk)<<"deltkk"<<deltkk;
            //std::cout<<"fq"<<fq(q,k,kl,kh)<<"et"<<et;
            //std::cout<<'u'<<random.u[i]<<"vv"<<random.v[i]<<std::endl;

        }
        //std::cout<<j<<"jj"<<x3<<"x3"<<b.u[j]<<"bu"<<b.v[j]<<"bv"<<std::endl;
        //std::cout<<1<<std::endl;
        //exit(1);
    }

    //std::ofstream fout;
    //fout.open ("matrixoutput.txt");

    double f0=edens(0.0);
	for(int j=0;j<nl;j++)
    {
        double x3=stepl*(j+0.5);
        double f=pow(edens(x3)/f0,yita);
        b.u[j]*=f;
        b.v[j]*=f;
        //fout<<x3<<' '<<b.u[j]<<' '<<b.v[j]<<' '<<sqrt(b.u[j]*b.u[j]+b.u[j]*b.u[j])<<endl;
        ////if (j > 10) break;
    }
    //这里的输出结果可以和Meyer environ.py 359附近B, psi = self._b.new_Bn()结果比较
    //fout.close();
    //exit(0);
    return b;
}


vec regb(double rr)
// Libanov et al. 2020
{
    if (rr>93) return 1.e-30;
    const double RR=93; //kpc
    const double tht=3.1415926/4; //45. deg
    //const double b0=8.3; //uG
    double cc=0.0599652;
//    using below code to calc cc
//    vec bb=regb(1.e-3);
//    double btt=sqrt(bb.v[0]*bb.v[0]+bb.v[1]*bb.v[1]+bb.v[2]*bb.v[2]);
//    cout<<8.3/btt<<endl;
    double r1=rr/RR;
    //al is the lowest non-0 root of tan(al)=3*al/(3-al*al);
    //namely pi/2.?
    double al=5.76; //see arXiv:1008.5353

    //double f0=cc*al*al*(al*cos(al)-sin(al));
    double f0=182.67*cc;
    double ff=cc*(al*cos(al*r1)-sin(al*r1)/r1)-f0*r1*r1/(al*al);
    double fp=cc*(-al*al*sin(al*r1)-al*cos(al*r1)/r1+sin(al*r1)/(r1*r1))-2*f0*r1/(al*al);
    double br=2*cos(tht)*ff/(r1*r1);
    double bth=-sin(tht)*fp/r1;
    double bph=al*sin(tht)*ff/r1;
    vec bb(br,bth,bph);
    //cout<<rr<<" "<<br<<" "<<bth<<" "<<bph<<endl;
    return bb;
}


field b_icm_reg()
//regular B field of Libanov2020
{
    field bb(nl,0.0); //set the length the same with turbulence B for better compatible with other parts of the codes
    std::ifstream ff;
    double x1[57],y1[57],x2[57],y2[57],x3[57],y3[57]; // 57 lines of data
    ff.open("Bicm_reg.txt");
    char szBuf[256];
    ff.getline(szBuf,sizeof(szBuf));
    for(int i=0; i<57; i++)
        ff>>x1[i]>>y1[i]>>x2[i]>>y2[i]>>x3[i]>>y3[i];
    //#r  Br   r   Btheta   r     Bphi
    vector<double> x1v(x1,x1+57);
    vector<double> x2v(x2,x2+57);
    vector<double> x3v(x3,x3+57);
    vector<double> y1v(y1,y1+57);
    vector<double> y2v(y2,y2+57);
    vector<double> y3v(y3,y3+57);

    for(int j=0; j<nl; j++)
    {
        double x3=stepl*(j+0.5); //每次都重新算一次x3其实很容易出错
        //bb.u[j]=interp1d(x2v,y2v,x3); //Btheta
        //bb.v[j]=interp1d(x3v,y3v,x3); //Bphi
        //cout<<x3<<" "<<bb.u[j]<<" "<<bb.v[j]<<endl;
        vec bb2(regb(x3));
        bb.u[j]=bb2.v[1];
        bb.v[j]=bb2.v[2];
        //cout<<j<<" "<<x3<<" "<<bb.u[j]<<" "<<bb.v[j]<<endl;
    }
    return bb; //uGauss
}



void faradayrm()
{
    arr rr=linspace(1.e-3,93,1000);
    arr res(1000);
    for(unsigned int j=0; j<rr.size(); j++)
    {
        vec bb = regb(rr[j]);
        res[j]=bb.v[0]*edens(rr[j]);
    }
    cout<<812.*quad(rr,res)<<endl; // regular B field of ICM
    //eq.(41) of https://ned.ipac.caltech.edu/level5/Sept04/Govoni/Govoni3_5.html

    const int nnn=1000;
    arr2 bb_field = set_bf_gcls(nnn,0); //其中一个垂向分量 nnn次realization
    arr rr2;
    arr res2(nl);
    arr rml(nl);
    for(int i=0; i<nl; i++) rr2.push_back(stepl*(i+0.5));
    double sum=0;
    for(int j=0; j<nnn; j++)
    {
        for(int i=0; i<nl; i++)
        {
            res2[i] = sqrt(2)*bb_field[j][i]*edens(rr2[i]);
            //eq.A8 of 1406.5972?
        }
        rml[j]=812.*quad(rr2,res2);
        cout<<rml[j]<<endl;
        sum += rml[j];
    }
    cout<<sum/nnn<<endl; // turbulent B field
}
