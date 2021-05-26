#include "alp.h"
#include "utils.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;


int test()
{
    double ma=30,ga=0.5;
    double ee = 1000000.; //1 GeV
    double bf = 1.; //
    double fi = 0.;
    double ne = 1.0;
    double ll = 100.; //kpc
    //cmat ddd = transmat(ma,ga,ee,bf,fi,ne,ll); //ll is coherent length
    //std::cout<<  <<std::endl;
    //ddd.print();
    transmat(ma,ga,ee,bf,fi,ne,ll).print();
    transmat(ma,ga,ee*10.,bf,fi,ne,ll).print();
    transmat(ma,ga,ee*100.,bf,fi,ne,ll).print();
    transmat(ma,ga,ee,bf*0.1,fi,ne,ll).print();
    transmat(ma,ga,ee,bf*100.,fi,ne,ll).print();
    transmat(ma,ga,ee,bf,3.1415926/4.,ne,ll).print();
    transmat(ma,ga,ee,bf,fi,ne,100000.).print();
    transmat(30.,0.5,10.1756333*1.e6,5.96377,3.92699,1.0,0.111408).print();
    transmat(30.,0.5,10.1756333*1.e6,5.96377,3.92699,0.043,0.111408).print();

    // fq(q,k,kl,kh)
    std::cout<<fq(-2.8,1.25521,0.1795,8.976)<<std::endl;
    std::cout<<fq(-2.8,0.1,0.1795,8.976)<<std::endl;
    std::cout<<fq(-2.8,5.0,0.1795,8.976)<<std::endl;
    std::cout<<fq(-2,0.1,0.1795,8.976)<<std::endl;
    std::cout<<fq(-2,5.0,0.1795,8.976)<<std::endl;
    std::cout<<fq(-3,0.1,0.1795,8.976)<<std::endl;
    std::cout<<fq(-3,5.0,0.1795,8.976)<<std::endl;

    return 1;
}

int par1()
{
    double gl=150.5793, gb=-13.2566, zred=0.017559;  // NGC 1275
    double ma=1,ga=0.3;
    char fi[]="rs_par1.txt";
    alp(ma,ga,gl,gb,zred,fi);
    return 1;
}

int scanmap()
{
    double gl=150.5793, gb=-13.2566, zred=0.017559;  // NGC 1275

    double ma0=0.1, ga0=0.01, ma1=100, ga1=100;
    double step_ma=(std::log10(ma1)-std::log10(ma0))/40.;
    double step_ga=(std::log10(ga1)-std::log10(ga0))/40.;

    for (int i=0; i<40; i++)
    {
        for (int j=0; j<40; j++)
        {
            double ma = pow(10.,log10(ma0)+step_ma*i);
            double ga = pow(10.,log10(ga0)+step_ga*j);

            char fi[100] = {};
            sprintf(fi, "ppp/rs_%d_%d.txt", i, j);
            cout<<i<<" "<<j<<" "<<ma<<" "<<ga<<" "<<endl;
            cout<<fi<<endl;
            alp(ma,ga,gl,gb,zred,fi);
        }
    }
    return 1;
}

int main()
{
    //test();
    //par1();
    //scanmap();
    faradayrm();
    return 1;
}


int alp(double ma, double ga, double gl, double gb, double zred, char* outfile)
{
	gl = gl*3.1415926/180.;
	gb = gb*3.1415926/180.;


    int nnn = 1; //5; //100;
    cout <<" l= "<< gl<<" b=" << gb << "  z="<<zred<< endl;
    cout<<"ma= "  << ma << "nev  ga=  " << ga << " 10^-11 GeV^-1， bnumber= " << nnn << "  " << endl;
    cout<<outfile<<endl;
    //std::exit(1);
    ////////////////////////////////////////////////////////

    vector<double> edata=getelist(); //GeV
    vector<double> tdata=readebl(edata,zred);
    //vector<double> tdata=readebl_dom11(edata,zred);

    //for(int i=0;i<edata.size();i++) cout<<edata[i]<<" "<<tdata[i]<<endl;
    //拷贝终端中的输出结果去文本中和python程序的ebl结果对比


    /////////////////////////////////////////////
    ////// 设置星系团磁场
//    vector<vector<double>> bb_field = set_bf_gcls(nnn,0);
//    vector<vector<double>> bb_field1 = set_bf_gcls(nnn,1); //bb_field, bb_field1分别是两个方向的分量
//    for(int k=0; k<nnn; k++)
//        for(int j=0; j<10; j++)
//            std::cout<<pow(bb_field[k][j]*bb_field[k][j]+bb_field1[k][j]*bb_field1[k][j],0.5)<<endl;

    //////////////////ICM规则磁场 用湍流时把下面的注释掉
    int regb=1;
    //if(regb){
    nnn = 1;
    field bb=b_icm_reg();
        vector<double> vec1(bb.u,bb.u+bb.size); //这里这样处理很别扭 后面把湍流、规则的输出结果类型统一起来
        vector<vector<double>> bb_field(1,vec1); //要么都vector要么都field
        vector<double> vec2(bb.v,bb.v+bb.size);
        vector<vector<double>> bb_field1(1,vec2);
    //}
//    for(int j=0; j<int(bb_field[0].size()); j++) {
//        double x3=stepl*(j+0.5);
//        cout<<x3<<" "<<bb_field[0][j]<<" "<<bb_field1[0][j]<<endl;
//    }
    /////////////ICM规则磁场


    /////////////////////////////////////////////
    ////// 计算传播
    int elen = edata.size();
    double *esave = new double[elen];
    double **psave = new double* [nnn];

    for (int k = 0; k < nnn; k++)
    {
        psave[k] = new double[elen]; //这样用大数组是不是很容易出现问题 用malloc或vector保险点 //不用 new是c++的 容易出错的是另一个c的
    }

    for (int k = 0; k < nnn; k++)
    {
        for (int j = 0; j < elen; j++)
        {
            if (j == 0 && (k % 24) == 0)
                cout << "k=" << k << "  j=" << 11 << endl;

            // 计算extragalactic background light 影响
            cmat transz(mat(1.0, 1.0, 1.0));
            transz.m[0][0] = sqrt(tdata[j]);
            transz.m[1][1] = sqrt(tdata[j]);
            cmat trans(mat(1.0, 1.0, 1.0));
            //edata[j] = 10.;
            trans = get_trans_gcls(k, edata[j], ma, ga, zred, bb_field, bb_field1);
            //cout<<"ma"<<ma<<"ga"<<ga<<"eee"<<edata[j]<<endl;
            //trans.print();

            cmat trans1(mat(1.0, 1.0, 1.0));
            trans1 = get_trans_gal(edata[j], ma, ga, gl, gb, zred);
            psave[k][j] = calc_pgg(trans1, transz, trans);
            esave[j] = edata[j];
        }
    }

    ////////////////////////////////////////////
    ////// 读出结果
    ofstream ofile;
    ofile.open(outfile);
    for(int j=0;j<elen;j++)
    {
        ofile<<esave[j]<<' ';
        for(int k=0;k<nnn;k++) ofile<<psave[k][j]<<' ';
        ofile<<endl;
    }
    ofile.close();

    ////  free the memory;
    for(int i=0;i<nnn;i++) delete[] psave[i];
    delete[] psave;
    delete[] esave;
    return 1;
}

//之前不知道怎么同时把bb_field1和bb_field返回出来，采用了一个笨办法return1
vector<vector<double>> set_bf_gcls(int nnn,int return1)
{
    //////////////////////////////////////////////////////////
    //随机数 for随机磁场
    std::mt19937 generator(103);
    //std::mt19937 generator (120);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    //////////////////////////////////////////////////////////

    double ffff;
    //double **bb_field = new double *[nnn];
    //double **bb_field1 = new double *[nnn];
    //for (int i = 0; i < nnn; i++)
    //{
    //    bb_field[i] = new double[nl];
    //    bb_field1[i] = new double[nl];
    //}
    vector<vector<double>> bb_field(nnn,vector<double>(nl));
    vector<vector<double>> bb_field1(nnn,vector<double>(nl));

    for (int k = 0; k < nnn; k++)
    {
        field random(nk, 0.0);
        field random1(nk, 0.0);
        cout<<k<<"-th realization. totally "<<nnn<<endl;

        for (int i = 0; i < nk; i++)
        {
            random.v[i] = dis(generator);
            ffff = dis(generator); //像是故意跳过一个随机数
            random.u[i] = dis(generator);
            ffff = dis(generator);
            random1.v[i] = dis(generator);
            ffff = dis(generator);
            random1.u[i] = dis(generator);
            ffff = dis(generator);
        }
        field b_source = b_fields(random, random1);
        for (int i = 0; i < nl; i++)
        {
            bb_field[k][i] = b_source.v[i]; //v u是磁场的两个分量 每个分量需要两个随机数
            bb_field1[k][i] = b_source.u[i];
        }
    }

    if (return1 == 1)
    {
        return bb_field1;
    } else {
        return bb_field;
    }

}

// 暂时先把nnn传进去 后面肯定要调整为nnn个模拟就调用bf_gcls函数nnn次就好
// int k---第k次模拟，目前这样处理主要是方便与原来的程序比较
cmat get_trans_gcls(int k, double eph, double ma, double ga, double zred,
        vector<vector<double>> &bb_field, vector<vector<double>> &bb_field1)
{
    /////////////////////////////////////////////
    ////// 计算星系团传播
    //
    double e = (1 + zred) * eph * 1.0e6;  //观测系到源系 // GeV 2 keV, why?--transmat函数的输入要求是keV, 见trans.cpp
    cmat trans(mat(1.0, 1.0, 1.0));
    for (int i = 0; i < nl; i++)
    {
        double r = stepl * (i+0.5);
        double b1 = bb_field[k][i];
        double b2 = bb_field1[k][i];
        double b = sqrt(b1 * b1 + b2 * b2);
//        double fi = atan(b2 / b1);
//        if (b1 < 0.0)
//            fi += pi;
        double fi=atan2(b1,b2);
        double ne = edens(r); ///////////////////////////
        //cout<<i<<" "<<e<<" "<<b<<" "<<fi<<" "<<ne<<" "<<stepl<<endl;
        cmat transtmp(transmat(ma, ga, e, b, fi, ne, stepl));
        trans = transtmp * trans;
        //transtmp.print();
        //cout<<" "<<endl;
	    //transmat(30.,0.5,10.1756333*1.e6,5.96377,3.92699,ne,0.111408).print();
	    // 可以与calcTransfer中self._Tn与return的值比较
    }
    /////////////////////////////////////////////
    //trans.print();
    //cout<<nl<<endl;
    return trans;
}

cmat get_trans_gal(double eph, double ma, double ga, double gl, double gb, double zred)
{
    //坐标；视线方向slight；垂直视线slight1，slight2
    vec slight = vec(cos(gb) * cos(gl), cos(gb) * sin(gl), sin(gb));
    vec earth(-8.5, 0, 0);
    vec ttp = slight + vec(1.0, 0.0, 0.0);
    vec slight1 = ttp - slight * (ttp * slight);
    slight1 = slight1 / slight1.mag();
    vec slight2 = slight.cross(slight1);

    /////////////////////////////////////////////
    double ee = (1 + zred) * eph * 1.0e6;  //观测系到源系 // GeV 2 keV, why?--transmat函数的输入要求是keV, 见trans.cpp
    ee = ee / (1 + zred);

    cmat trans1(mat(1.0, 1.0, 1.0));  //对角矩阵
    /////////////////////////////////////////////
    ////// 计算银河系传播
    //transform in Milk Way;
    for (int i = 0; i < ns; i++) //ns=100
    {
        vec posi = earth + slight * (i * mashm);
        vec b = b_fieldm(posi);
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
        trans1 = transmat(ma, ga, ee, bt, fi, ne, mashm) * trans1;
    }
    return trans1;
}


double calc_pgg(cmat &trans1, cmat &transz, cmat &trans)
{
    ////// 收集结果
    //trans1--银河系传播  transz--EBL吸收  trans--星系团传播
    cmat transtt = trans1 * transz * trans;
    //trans.print();
    //trans1.print();
    //cmat p0(mat(0, 0, 1)), p1(mat(0, 0, 1)); // for pro_aa
    //cmat pro1 = p1 * transtt * p0 * transtt.ct(); // for pro_aa
    cmat p0(mat(0.5, 0.5, 0)), p1(mat(1.0, 0.0, 0.0)), p2(mat(0.0, 1.0, 0.0));
    cmat pro1 = p1 * transtt * p0 * transtt.ct();
    cmat pro2 = p2 * transtt * p0 * transtt.ct();
    double progg = pro1.tr().re + pro2.tr().re;  //re是取实部
    //double progg = pro1.tr().re; // pro_aa
    return progg;
}


vector<double> getelist()
{
    int elen = eol;
    double e1=eo1, e2=eo2; //GeV
    double estep = (std::log10(e2)-std::log10(e1))/(elen-1);
    vector<double> edata(elen);

    for (int i = 0; i < elen; i++)
        edata[i] = pow(10.,std::log10(e1)+(i)*estep);
    return edata;
}


