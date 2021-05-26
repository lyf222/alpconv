#include "alp.h"
#include <fstream>
#include <cstring>
using namespace std;

int main(int argc, char **argv)
{
    double smaxm = 30;
    double gl=0, gb=0;
    double ma = 3.6, ga = 1.01;
    int bfiledFlag = 1;

    gl = std::atof(argv[1]);
    gb = std::atof(argv[2]);
    smaxm = std::atof(argv[3]);
    ma = std::atof(argv[4]);
    ga = std::atof(argv[5]);
    bfiledFlag = std::atoi(argv[6]);
    //////////////////////////////////////////////////////////
    cout <<" l= "<< gl<<" b=" << gb << "  dis="<<smaxm<< endl;
    cout <<"ma = "<< ma<<" nev" << "      ga = " << ga<<" x 10^-11 GeV ^-1" << endl;
    cout << "bfiledFlag:  " << bfiledFlag << endl;
    
    //alp(smaxm,gl,gb,ma,ga,bfiledFlag);
    return 1;
}

int alp(double smaxm, double gl, double gb, double ma, double ga, int bfiledFlag, 
		double *energies, int nenergies, double *prob_out)
{
    gb=gb/180.0*pi;
    gl=gl/180.0*pi;
    double mashm = smaxm / ns;
    vec slight = vec(cos(gb) * cos(gl), cos(gb) * sin(gl), sin(gb));
    vec earth(-8.5, 0, 0);
    vec ttp = slight + vec(1.0, 0.0, 0.0);
    vec slight1 = ttp - slight * (ttp * slight);
    slight1 = slight1 / slight1.mag();
    vec slight2 = slight.cross(slight1);
    ////////////////////////////////////////////////////////
    int elen = nenergies;
    double *edata = new double[elen];
    for (int i = 0; i < nenergies; i++)
    {
        edata[i] = energies[i];
    }
    ///////////////////////////////////////////////
    /////////////////////////////////////////////
    double esave[elen], psave[elen];
    for (int j = 0; j < elen; j++)
    {
        double e = edata[j] * 1.0e6;

        cmat trans1(mat(1.0, 1.0, 1.0));
        //transform in Milk Way;
        for (int i = 0; i < ns; i++)
        {
            vec posi = earth + slight * (i * mashm);
	    vec b(0);
	    if (bfiledFlag == 1) b = b_fieldm1(posi);
	    if (bfiledFlag == 2) b = b_fieldm2(posi);
	    if (bfiledFlag == 3) b = b_fieldm(posi);
            double bt1 = b * slight1;
            double bt2 = b * slight2;
            double bt = sqrt(bt1 * bt1 + bt2 * bt2);
            double fi = atan(bt2 / bt1);
            if (bt1 < 0.0)
                fi += pi;
            if (bt1 == 0)
                fi = pi;
            if (bt == 0)
                continue;
            double ne = edenmw(posi.v[0], posi.v[1], posi.v[2]);
            trans1 = transmat(ma, ga, e, bt, fi, ne, mashm) * trans1;
        }

        cmat transtt = trans1;
        cmat p0(mat(0.5, 0.5, 0)), p1(mat(1.0, 0.0, 0.0)), p2(mat(0.0, 1.0, 0.0));
        cmat pro1 = p1 * transtt * p0 * transtt.ct();
        cmat pro2 = p2 * transtt * p0 * transtt.ct();
        double progg = pro1.tr().re + pro2.tr().re;
        esave[j] = e * 1.0e-6;
        {
            psave[j] = progg;
        }
    }
    for (int j = 0; j < elen; j++)
    {
        prob_out[j] = psave[j];
    }
    return 1;
}
