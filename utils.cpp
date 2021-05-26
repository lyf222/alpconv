#include "utils.h"

using namespace std;

void printv(const vector<double>& v)
{
    for(vector<double>::const_iterator it = v.begin(); it!=v.end(); it++)
        cout<<*it<<" ";
    cout << endl;
}

vector<double> linspace(double start, double stop, int num)
//已对比 与np.linspace相同
{
    vector<double> v(num);
    double bins=(stop-start)/(num-1);
    for(int i=0; i<num; i++) v[i]=start+i*bins;
    return v;
}

vector<double> arange(double start, double stop, double step)
//已对比 与np.arange相同
{
    int num=(stop-start)/step;
    vector<double> v(num);
    for(int i=0; i<num; i++) v[i]=start+step*i;
    return v;
}

void split(const string& s, vector<string>& tokens, const string& delimiters)
//等价于python中的s.split()
//默认参数在函数声明中提供，当又有声明又有定义时，定义中不允许默认参数
{
    string::size_type lastPos = s.find_first_not_of(delimiters, 0);
    string::size_type pos =s.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        tokens.push_back(s.substr(lastPos, pos - lastPos));
        lastPos = s.find_first_not_of(delimiters, pos);
        pos = s.find_first_of(delimiters, lastPos);
    }
}


double quad(vector<double> &xx, vector<double> &yy)
//simplest trapezoid integration with given sample points
//通过自己多运行一遍quad函数 将采样点翻倍来判断精度
{
    double sum=0;
    for(unsigned int i=1;i<yy.size();i++) {
        sum+=(yy.at(i)+yy.at(i-1))*(xx.at(i)-xx.at(i-1))/2.;
    }
    return sum;
}
//vector<double> xx=linspace(0.,4.,100);
//vector<double> yy(xx.size());
//for(int i=0;i<yy.size();i++) yy[i]=pow(xx[i],2);
//cout<<xx.size()<<" "<<quad(xx,yy)<<endl;
//xx=linspace(0.,4.,200);
//yy.resize(xx.size());
//for(int i=0;i<yy.size();i++) yy[i]=pow(xx[i],2);
//cout<<xx.size()<<" "<<quad(xx,yy)<<endl;

double interp1d(vector<double> &xvec, vector<double> &yvec, double xx)
//如果不做外插 只是线性内插 可直接调用该函数
//作为插值函数的模板 以后要用插值 拷贝过去在此基础上改即可
{
    if(xx<=xvec[0]) return yvec[0]; //改这里决定往前外插的值
    int xsi=xvec.size();
    for(int i=1; i<xsi; i++) //插值
    {
        if(xvec[i]>xx)
            return yvec[i-1]+(yvec[i]-yvec[i-1])/(xvec[i]-xvec[i-1])*(xx-xvec[i-1]);
    }
    return 1.e-30; //yvec.back(); //改这里决定往后外插的值
}

vector<vector<double>> readtxt(const char* filename, int ncols, int nSkip)
{
    vector<double> vec;
    vector<vector<double>> data;

    ifstream ff;
    ff.open(filename);
    //skip first nSkip lines
    char szBuf[512];
    for(int i=0; i<nSkip; i++) ff.getline(szBuf,sizeof(szBuf));
    //read all data
    double tmp;
    while(ff>>tmp) vec.push_back(tmp);
    ff.close();
    //reshape the data
    int nrows=vec.size()/ncols;
    data.resize(ncols);
    for(int j=0; j<ncols; j++) data[j].resize(nrows);
    for(int j=0; j<ncols; j++)
        for(int i=0; i<nrows; i++)
            data[j][i]=vec.at(i*ncols+j);
    //for(int j=0;j<nrows;j++) cout<<data[0][j]<<endl; //print first line
    return data;
}
