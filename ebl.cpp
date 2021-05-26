#include "alp.h"
#include "utils.h"

using namespace std;

vector<double> readebl(vector<double> &edata,double redshift)
{
    ///////////////////////////////////////
    // 读extragalactic background light
    vector<double> tdata(edata.size());
    double xt[50], yt[50];
    field dd = readebl_fran08(redshift);

    for(int i=0; i<50; i++)
    {
        xt[i]=dd.u[i];
        yt[i]=dd.v[i];
    }

    for(unsigned int i=0; i<edata.size(); i++)
    {
        double e=edata[i]/1000;
        if(e<xt[0])
        {
            tdata[i] = 1.0;
            continue;
        }
        for (int j=1; j<50; j++) //插值
        {
            if (xt[j]>e)
            {
                tdata[i] = yt[j-1] + (yt[j] - yt[j-1]) / (xt[j] - xt[j-1]) * (e - xt[j-1]);
                tdata[i] = 1.0 / tdata[i]; //e^tau 2 e^-tau
                break;
            }
        }
    }
    return tdata;
}


field readebl_fran08(double zz)
//return res.u: E(TeV); res.v: e^tau
{
    if(zz<1.e-3) std::exit(0); //checking whether the zz is within the range

    string temp;
    double tmpp,zt1=0,zlas=0;
    field res(50,0),las(50,0);
    ifstream xytfile;

    xytfile.open("ebl/Gamma-gamma-opacity-z=0-1.dat");
    //https://arxiv.org/pdf/0805.1841.pdf
    //http://www.astro.unipd.it/background/
    //http://www.astro.unipd.it/background/tev/Gamma-gamma-opacity-z=0-1.dat.gz

    for(int i=0;i<1000;i++) //from 1.e-3 to 1.001000, one z point per dz=0.001
    {
        for(int j=0;j<50;j++){ //store the last one
            las.u[j]=res.u[j];
            las.v[j]=res.v[j];
        }
        zlas = zt1;

        vector<string> tokens;
        getline(xytfile,temp);
        split(temp,tokens);
        zt1 = atof(tokens[3].c_str());
        getline(xytfile,temp);
        // cout<<zt1<<endl;

        for(int j=0; j<50; j++)
        {
            getline(xytfile,temp);
            istringstream stream(temp);
            stream>>res.u[j]>>tmpp>>tmpp>>res.v[j];
            //E0 [TeV],    E0 [eV],   tau(E0)   e^tau
        }
        if(abs(zt1-zz)<0.00001) break;  //return res;
        if(zt1>zz)   //找bin的右边界
        {
            //zt1,res.v; //当前
            //zlas,las.v; //前一个
            //y = (x - x1)/(x2 - x1)*(y2 - y1) + y1;
            for(int j=0;j<50;j++)
                res.v[j]=las.v[j]+(res.v[j]-las.v[j])*((zz-zlas)/(zt1-zlas));
            break; //return res;
        }
    }
    return res;
}

vector<double> readebl_dom11(vector<double> &edata,double redshift)
// aims to for a list of energies and a given redshift
// return corresponding e^-tau
{
    vector<double> zl= {0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895,
                        0.11684211, 0.13210526,  0.14736842,    0.16263158, 0.17789474,  0.19315789,
                        0.20842105,  0.22368421,  0.23894737, 0.25421053,  0.26947368,  0.28473684,
                        0.3,  0.35,  0.4, 0.45,  0.5,  0.55,  0.6,  0.65,  0.7,  0.75,  0.8,  0.85,
                        0.9,   0.95,  1.,1.2,1.4,1.6,1.8,2.
                       };
    if(redshift<zl[0])
    {
        cout<<"error!";    //checking whether the zz is within the range
        std::exit(0);
    }

    vector<vector<double>> data=readtxt("ebl/tau_dominguez11.out",40,5);
    vector<double> res(data[0].size(),0);

    //下面的内容可以提取成一个2d插值函数
    for(unsigned int j=1; j<zl.size(); j++)
    {
        //zl[j],data[j+1] //now //for data its 1st col is energy
        //zl[j-1],data[j] //last
        if(zl[j]>redshift)
        {
            for(unsigned int i=0; i<data[0].size(); i++)
                res[i]=data[j][i]+(data[j+1][i]-data[j][i])*((redshift-zl[j-1])/(zl[j]-zl[j-1])); //tau
            break;
        }
    }
    for(unsigned int i=0; i<data[0].size(); i++) cout<<data[0][i]<<" "<<res[i]<<endl;
    //for sorted edata, the way used in readebl function is faster
    vector<double> etau(edata.size(),0);
    for(unsigned int i=0; i<edata.size(); i++)
    {
        etau[i]=exp(-interp1d(data[0],res,edata[i]/1.e3));
        //cout<<i<<" "<<edata[i]<<" "<<etau[i]<<" "<<edata.size()<<endl;
    }

    return etau;
}

//using c++ to rewrite code of dealing with the ebl data is bothering.
//since we do not consider the IGM Bfield
//it is good approximation of pgg=p_icm,gg*ebl*p_mw,gg+p_icm,ga*p_mw,ag
//that we can use the python ebl code
