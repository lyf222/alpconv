#include"alp.h"
#include "math.h"
using namespace std;
//electron density in Milky Way and Gamma-ray source
inline double sech2(double x)
{
	return (1.0 / cosh(x))*(1.0 / cosh(x));
}
double edenmw(double x,double y, double z)
{
	double n1=0.033/0.95, h1=0.95, a1=17;
	double n2=0.09,h2=0.14, a2=3.7;
	double r = sqrt(x*x + y*y);
	double g1;
	if(r>a1)
    {
        g1=0.0;
    }
    else
    {
        g1=cos(pi/2*r*a1)/cos(pi/2*rsun*a1);
    }
    double tmp=(r-a2)/a2;
    double g2=exp(-tmp*tmp);
	return n1*g1*sech2(z/h1)+n2*g2*sech2(z/h2);
}
double edens(double r)
//r: kpc?  return: cm^-3?
{
    return 3.9e-2/(pow((1+r*r/80/80),1.8))+4.05e-3/(pow((1+r*r/280/280),0.87));
}
