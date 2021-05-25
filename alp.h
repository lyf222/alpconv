#ifndef alp_h
#define alp_h
#include<vector>
#include<iostream>
#include<math.h>
#include<random>
#include<time.h>
#include <sstream>
#include <fstream>

const double pi=3.1415926;
const double rsun=8.5;
// mash the sight line in Milk  Way
const int ns=100;
const double smaxm=30;
const double mashm=smaxm/ns;
//mash the engergy
const int nodee=100;
const double masse=0.511e-3; //electron mass in Gev

//parameter  for Milky Way magnetic field
const double sdisk[8]={5.1,6.3,7.1,8.3,9.8,11.4,12.7,15.5};
//const double bd[8]={0.1,3.0,-0.9,-0.8,-2.0,-4.2,0.0,2.7}; //original parameters
const double bd[8]={0.1,3.0,-0.9,-0.8,-2.0,-3.5,0.0,2.7}; //planck parameters, bd6=-4.2 --> bd6=-3.5
const double openi=11.5/180*pi,bring=0.1;

const double hdisk=0.4,wdisk=0.27;
const double bn=1.4,bs=-1.1,rn=9.22,rs=16.7,wh=0.20,z0=5.3;
//const double bx=4.6,theta0=49.0/180*pi,rxc=4.8,rx=2.9; //orignial parameters
const double bx=1.8,theta0=49.0/180*pi,rxc=4.8,rx=2.9; //planck parameters, bx=4.6 --> bx=1.8

//parameter for magnetic field in 星系团磁场
//const double kl=2*pi/35,kh=2*pi/0.7,q=-2.8,rmax=370,sigamb=10.0,yita=0.5;
//const double kl=2*pi/35,kh=2*pi/0.7,q=-2.8,rmax=500.,sigamb=10.0,yita=0.5;
const double kl=0.1795,kh=8.976,q=-2.8,rmax=500.,sigamb=10.0,yita=0.5;


//const int nk=int(20*(log10(kh/kl)+3));  // nk is the total number of spacings in k; where k is the wave number
//const int nl= (rmax*kh)/4; // l is the integration step length, thus nl is the number of step

const int nl= (rmax*kh)/1; // l is the integration step length, thus nl is the number of step
const double lmin=kl*1.0e-3,stepl=rmax/nl;
const int nk= 10.*(log10(kh)-log10(lmin))*(log10(kh)-log10(lmin))+1; //for comparing with Meyer's code
const double stepk=pow(10,(log10(kh/kl)+3)/nk);

//const double eo1=10,eo2=pow(10,3.5); //GeV
const double eo1=0.1,eo2=1000; //GeV
const int eol = 250; //used in getelist()





inline void Qsort(double* a, int low, int high)
{
    if(low >= high)
    {
        return;
    }
    int first = low;
    int last = high;
    double key = a[first];/*用字表的第一个记录作为枢轴*/

    while(first < last)
    {
        while(first < last && a[last] >= key)
        {
            --last;
        }

        a[first] = a[last];/*将比第一个小的移到低端*/

        while(first < last && a[first] <= key)
        {
            ++first;
        }

        a[last] = a[first];
/*将比第一个大的移到高端*/
    }
    a[first] = key;/*枢轴记录到位*/
    Qsort(a, low, first-1);
    Qsort(a, first+1, high);
}
class field
{
    public:
    int size;
    double *u;
    double *v;
    field()
    {
        size=0;
        u=NULL;
        v=NULL;
    }
    field(int s,double c)
    {
        size=s;
        u=new double[s];
        v=new double[s];
        for(int i=0;i<s;i++)
        {
            u[i]=c;
            v[i]=c;
        }
    }
    ~ field()
    {
        delete[] u;
        delete[] v;
    }

};

class field3
{
    public:
    int size;
    double *u;
    double *v;
    double *w;
    field3()
    {
        size=0;
        u=NULL;
        v=NULL;
        w=NULL;
    }
    field3(int s,double c)
    {
        size=s;
        u=new double[s];
        v=new double[s];
        w=new double[s];
        for(int i=0;i<s;i++)
        {
            u[i]=c;
            v[i]=c;
            w[i]=c;
        }
    }
    ~ field3()
    {
        delete[] u;
        delete[] v;
        delete[] w;
    }

};

class vec
{
public:
    double v[3];
    vec(double x, double y, double z)
    {
        v[0]=x;
        v[1]=y;
        v[2]=z;
    }
    vec(double c)
    {
        v[0]=c;
        v[1]=c;
        v[2]=c;
    }
    vec operator*(double b)
    {
        vec c(0.0);
        for(int i=0;i<3;i++) c.v[i]=v[i]*b;
        return c;
    }
    vec operator/(double b)
    {
        vec c(0.0);
        for(int i=0;i<3;i++) c.v[i]=v[i]/b;
        return c;
    }
    double operator*(vec b)
    {
        double c=0;
        for(int i=0;i<3;i++) c+=v[i]*b.v[i];
        return c;
    }
    vec  operator+(vec b)
    {
        vec c(0.0);
        for(int i=0;i<3;i++) c.v[i]=v[i]+b.v[i];
        return c;
    }
    vec  operator-(vec b)
    {
        vec c(0.0);
        for(int i=0;i<3;i++) c.v[i]=v[i]-b.v[i];
        return c;
    }
    vec cross(vec b)
    {
        vec c(0.0);
        c.v[0]=v[1]*b.v[2]-v[2]*b.v[1];
        c.v[1]=v[2]*b.v[0]-v[0]*b.v[2];
        c.v[2]=v[0]*b.v[1]-v[1]*b.v[0];
        return c;
    }
    double mag()
    {
        return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    }
    void print()
    {
        std::cout<<v[0]<<"  "<<v[1]<<"  "<<v[2]<<std::endl;
    }

};
class mat
{
public:
    double m[3][3];
    mat(double m11,double m12,double m13,double m21,double m22,double m23,double m31,double m32,double m33)
    {
        m[0][0]=m11;
        m[0][1]=m12;
        m[0][2]=m13;
        m[1][0]=m21;
        m[1][1]=m22;
        m[1][2]=m23;
        m[2][0]=m31;
        m[2][1]=m32;
        m[2][2]=m33;
    }
    mat(double m11,double m22,double m33)
    {
        m[0][0]=m11;
        m[0][1]=0.0;
        m[0][2]=0.0;
        m[1][0]=0.0;
        m[1][1]=m22;
        m[1][2]=0.0;
        m[2][0]=0.0;
        m[2][1]=0.0;
        m[2][2]=m33;
    }
    mat(double* mm)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            m[i][j]=mm[i*3+j];
    }
    mat(double c)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            m[i][j]=c;
    }
    mat operator*(mat b)
    {
        mat c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                    c.m[i][j]+=m[i][k]*b.m[k][j];
        return c;

    }
    mat operator*(double b)
    {
        mat c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[i][j]*b;
        return c;

    }

    vec operator*(vec b)
    {
        vec c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.v[i]+=m[i][j]*b.v[j];
        return c;
    }
    mat t()
    {
        mat c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[j][i];
        return c;
    }
    void print()
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                std::cout<<m[i][j]<<"  ";
            }
            std::cout<<std::endl;
        }
    }

};
class cplex
{
public:
    double re;
    double im;
    cplex(){};
    cplex(double a):re(a),im(0.0){};
    cplex(double a,double b):re(a),im(b){};
    cplex operator*(cplex b)
    {
        return cplex(re*b.re-im*b.im,re*b.im+im*b.re);
    }
    cplex operator*(double b)
    {
        return cplex(re*b,im*b);
    }
    cplex operator+(cplex b)
    {
        return cplex(re+b.re,im+b.im);
    }

    cplex operator-(cplex b)
    {
        return cplex(re-b.re,im-b.im);
    }
    void print()
    {
        std::cout<<"  "<<re<<"  "<<im<<"i   ";
    }
    cplex conj()
    {
        return cplex(re,-im);
    }
};

inline cplex cexp(cplex a)
{
    cplex b;
    b.re=exp(a.re)*cos(a.im);
    b.im=exp(a.re)*sin(a.im);
    return b;
}
class cmat
{
    public:
    cplex m[3][3];
    cmat();
    cmat(cplex m11,cplex m12,cplex m13,cplex m21,cplex m22,cplex m23,cplex m31,cplex m32,cplex m33)
    {
        m[0][0]=m11;
        m[0][1]=m12;
        m[0][2]=m13;
        m[1][0]=m21;
        m[1][1]=m22;
        m[1][2]=m23;
        m[2][0]=m31;
        m[2][1]=m32;
        m[2][2]=m33;
    }
    cmat(cplex m11,cplex m22,cplex m33)
    {
        m[0][0]=m11;
        m[0][1]=cplex(0.0);
        m[0][2]=cplex(0.0);
        m[1][0]=cplex(0.0);
        m[1][1]=m22;
        m[1][2]=cplex(0.0);
        m[2][0]=cplex(0.0);
        m[2][1]=cplex(0.0);
        m[2][2]=m33;
    }
    cmat(cplex* mm)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            m[i][j]=mm[i*3+j];
    }
    cmat(mat c)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            m[i][j]=cplex(c.m[i][j]);
    }
    cmat(cplex c)
    {
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
            m[i][j]=c;
    }
    cmat operator*(cmat b)
    {
        cmat c(cplex(0.0));
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                    c.m[i][j]=c.m[i][j]+m[i][k]*b.m[k][j];
        return c;

    }
    cmat operator+(cmat b)
    {
        cmat c(cplex(0.0));
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[i][j]+b.m[i][j];
        return c;

    }
    cmat operator*(cplex b)
    {
        cmat c(cplex(0.0));
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[i][j]*b;
        return c;
    }
    friend cmat operator*(mat a,cmat b)
    {
        cmat c(cplex(0.0));
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                for(int k=0;k<3;k++)
                    c.m[i][j]=c.m[i][j]+cplex(a.m[i][k])*b.m[k][j];
        return c;
    }
    void print()
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                m[i][j].print();
            }
            std::cout<<std::endl;
        }
    }
    mat get_re()
    {
        mat c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[i][j].re;
        return c;
    }
    mat get_im()
    {
        mat c(0.0);
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[i][j].im;
        return c;
    }
    cmat ct()
    {
        cmat c(cplex(0.0));
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                c.m[i][j]=m[j][i].conj();
        return c;
    }
    cplex tr()
    {
        return m[0][0]+m[1][1]+m[2][2];
    }

};


using arr=std::vector<double>;
using arr2=std::vector<std::vector<double>>;

int alp(double ma, double ga, double gl, double gb, double zred, char* outfile);
double edenmw(double x,double y, double z);
double edens(double r);
vec b_fieldm(vec cordxy);
field b_fields(field& random,field& random1);
cmat transmat(double ma,double ga,double e,double b,double fi,double ne,double l);
vec cord(double ra,double dec);
double fq(double q,double k,double kl,double kh);
double fq1(double q,double s1,double s2);
double integ(double (*func)(double c),double low,double high,int node=100);
double integnoise(double (*func)(double c,double n),double low,double high,double nose,int node=100);
double integnoise2(double (*func)(double c,double n,double n1),double low,double high,double noise,double noise1,int node=100);
double interp(double *x,double *y,int node,double xx);
vec b_fieldm2(vec cordxy);
vec b_fieldm1(vec cordxy);
vec b_fieldm3(vec cordxy);
field b_icm_reg();
vec regb(double rr);
void faradayrm();


double calc_pgg(cmat &trans1, cmat &transz, cmat &trans);
arr2 set_bf_gcls(int nnn,int return1);
arr getelist();
arr readebl(arr &edata,double redshift);
field readebl_fran08(double zz);
arr readebl_dom11(arr &edata,double redshift);

cmat get_trans_gcls(int k, double eph, double ma, double ga, double zred,
        arr2 &bb_field, arr2 &bb_field1);
cmat get_trans_gal(double eph, double ma, double ga, double gl, double gb, double zred);

#endif
